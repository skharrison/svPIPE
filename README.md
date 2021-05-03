# svPIPE

### Structural Variant Analysis 

**Programs used for SV calling:**
- [Pindel](https://github.com/genome/pindel) Helpful documentation [here](http://gmt.genome.wustl.edu/packages/pindel/user-manual.html)
- [Breakdancer](https://github.com/genome/breakdancer) Helpful documentation [here](https://gmt.genome.wustl.edu/packages/breakdancer/documentation.html) 
- [DELLY](https://github.com/dellytools/delly)
- [Manta](https://github.com/Illumina/manta)
- [LUMPY](https://github.com/arq5x/lumpy-sv)
- [Gridss](https://github.com/PapenfussLab/gridss)

**Other Dependencies**
- [Samtools](https://github.com/samtools/samtools)
- [Bcftools](https://github.com/samtools/bcftools)
- [Bedtools](https://github.com/arq5x/bedtools2)
- [IGV](https://github.com/igvteam/igv)
- [Jupyter Notebook](https://github.com/jupyter/notebook) (optional but example pipeline inside of notebook but requires SOS notebook from jupyterlab to use)
- [Stuctural Variant Annotation](https://www.bioconductor.org/packages/release/bioc/html/StructuralVariantAnnotation.html) R package used to annotate Gridss and Lumpy output

#### How to run pipeline after software installation 

- Clone repo and add scripts to data folder 
- run the runScript.sh like so:
```
svFolder= path/to/sv/foldertestSV            
cd $svFolder            
chmod u+x runScript.sh
mydata=$(pwd)
refs=$mydata/refs           ##replace refs with reference folder name 
bams=$mydata/bams           ##replace bams with bam folder name
reference=$refs/refSeqName  ##replace refSeqName with name of reference file inside of refs folder 
insertSize = 250            ##input estimated insert size based on sequence data 

./runScript.sh -d $mydata -b $bams -r $reference -i 250
```
- convert delly output using bcftools:
```
##converting delly from bcf to vcf 

cd $svFolder/delly
for file in *.bcf ;
do
out=${file%.bcf}.vcf
bcftools view -Ov -o $out $file
done
```
- Run reformat python script 
```
for b in $svFolder/stats/*.stats
do
tmp=${b##*/}
sample=${tmp%.sorted.stats}
python3 $svFolder/reformat_merge.py $projectD/$svFolder $sample
done
```
- Overlap all samples calls after reformatted
```
numSamples=20                   ## replace 20 with number of samples 
python3 $projectD/$svFolder/overlapSamples.py $projectD/$svFolder/mergeSV $numSamples
```

### Coverage and verification of merged regions:

#### Structural Variant Analysis GUI

* SVGui.java is a graphic interface to batch view igv images and automate other portions of pipeline (***still a work in progress***)
* Quick youtube video to demo what GUI looks like and what each tool will display can be found [here](https://www.youtube.com/watch?v=kPWZuFNhOJI&feature=youtu.be)

#### IGV Displayer Tool

##### HOW TO USE:
---------------
  * First upload image headers using **Add Strain Labels** button. Images labels should be in the same order as they appear right to left in IGV image.  
  * Then upload all images (must be .png,.jpg, or, .svg) images uing **Browse** button 
  * Once all images of regions wish to display inside of IGV are loaded click **View Images** There will be a small popup indicating that image rendering and         table building is occuring. When this is finished the table will appear and can click out of the popup. 
  * Can now batch view all the regions where the CHROM, START, STOP, are on the left side of the table and the image is on the right. If an image looks as if it       is of interest can then click the checkbox. After looking through the photos can then choose to save all of these regions that appear to be of interest by         hitting the **Save Checked** button. Generating a new bed file for further analysis if desired. 
  
##### NOTES ABOUT IMAGE FORMAT AND GENERATION:
----
Images were generated of desired regions by first taking all regions (in bed format) in which would like to validate/analyze and running bedtools like so:
- Make a folder where you would like images to be, in example below made folder igvFolder
- Makes bash script to run inside IGV for image generation:
```
bedtools igv -path path/to/igvFolder -slop 100 -i $svFolder/mergeSV/noShareSV.bed > $svFolder/allSV_igvScript.sh
```
> slop is the number of bp would like to be flanking region (must include for proper import into GUI so if wish to not have any just put 0)

* All sample bams then loaded into IGV (set preferences in IGV to COLLAPSED view) can move samples around to order would like them to be in GUI table (top to bottom) will be equivalent to right to left.
* Usually go to view and remove Name and Header panels to produce a cleaner looking image to display inside table 
* Click igvTools ---> run batch script  (to automatically generate images of all regions
***NOTE:*** Changed alignment preferences to be able to visualize larger regions then the default of 20kp. However with a large number of samples and visualizing large regions is dependent on  memory capacity. 
* Once all images have been produced by igvTools then flip horizontally. Easy to batch do in a folder by using [magick](https://imagemagick.org/script/mogrify.php)
```
for photo in *.png ; do convert $photo -rotate 90 $photo ; done
```
* Images should now be ready to load into GUI

#### Compute Coverage
* Tool to run coverage analysis of a set of sample bams for a specific bam file. 
* Purpose is to be able to determine coverage of regions that might have SV in order as this would likely influence coverage in those regions for samples           containing the variant.
* This tool is running the bedtools multicov command. (**MUST HAVE BEDTOOLS PRE-INSTALLED**). 

##### HOW TO USE:
-----------
* Load in bams of interest (*Note: import order of files will be order of samples in output file).
* Load in the regions that would like to compute coverage in this bams (in bed format). 
* Choose output file for the coverage calculations to be saved to. 
* Once all necessary files are imported can then click Submit
*bar will appear showing that bedtools is processing. When bedtools has finished the bar will turn to done and output file will be displayed on the screen. 

***NOTE: BEDTOOLS is not compatible with Windows so this tool will NOT work on Windows machines*** 

#### Color Sample Table 

* Input a table containing the variants, all samples, and normalized coverage values. Given a choice of highlighting color (with default color yellow), the tool highlights samples that have the specific variant.
* Input table generated by structural variant pipeline mentioned above. 
