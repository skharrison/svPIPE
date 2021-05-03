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

#### Structural Variant Analysis GUI

* SVGui.java is a graphic interface to batch view igv images and automate other portions of pipeline (***still a work in progress***)
* Quick youtube video to demo what GUI looks like and what each tool will display can be found [here](https://www.youtube.com/watch?v=kPWZuFNhOJI&feature=youtu.be)
* JAVA GUI code and information available at: https://github.com/skharrison/progFinal


##TODO:
- SET UP SOS jupyter notebook to easily run the entire pipeline 
