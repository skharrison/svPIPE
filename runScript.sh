#!/bin/bash


while getopts d:b:r:i: flag
do
    case "${flag}" in
        d) outDir=${OPTARG};;
        b) bams=${OPTARG};;
        r) reference=${OPTARG};;
        i) insert=${OPTARG};;
    esac
done
echo "Data located in: $outDir"
echo "Bams files located in: $bams";
echo "Reference file: $reference";
echo "Insert size: $insert";


##moves to data output folder and makes out directories 
cd $outDir
mkdir -p pindel
mkdir -p delly
mkdir -p breakdancer
mkdir -p manta
mkdir -p gridss
mkdir -p lumpy 
mkdir -p stats



for b in $bams/*.bam
do

tmp=${b##*/}
stats=${tmp%.bam}.stats
#echo "samtools stats $b > $outDir/stats/$stats"
samtools stats $b > $outDir/stats/$stats

#running breakdancer
breakCfg=${tmp%.bam}.break.cfg
breakOut=${tmp%.sorted.bam}_breakdancer.txt
perl /usr/local/breakdancer/perl/bam2cfg.pl -g -h $b > $outDir/breakdancer/$breakCfg
breakdancer-max -r 4 $outDir/breakdancer/$breakCfg > $outDir/breakdancer/$breakOut

##setup delly 
dellyOut1=${tmp%.sorted.bam}.bcf
delly call -g $reference -o $outDir/delly/$dellyOut1 $b

##running manta  
sample=${tmp%.sorted.bam}
mkdir $outDir/manta/$sample
configManta.py --bam $b --referenceFasta $reference --runDir $outDir/manta/$sample
$outDir/manta/$sample/runWorkflow.py

##running gridss
gridOut=${tmp%.sorted.bam}_gridss.vcf.gz
gridAs=${tmp%.sorted.bam}_g.bam
gridss.sh --jar /usr/local/gridss/scripts/gridss.jar --reference $reference --output $outDir/gridss/$gridOut --labels $sample --assembly $outDir/gridss/$gridAs $b

##setup lumpy
discord=${b%.sorted.bam}.discordants.bam
splitters=${b%.sorted.bam}.splitters.bam
histo=${b%.sorted.bam}.lib1.histo

samtools view -b -F 1294 $b > $discord
samtools view -h $b \
    | /usr/local/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > $splitters
samtools view $b \
    | tail -n+100000 \
    | /usr/local/lumpy-sv/scripts/pairend_distro.py \
    ###TODO
    -r $insert \
    -X 4 \
    -N 10000 \
    -o $histo
done

##run lumpy in background
python3 $outDir/makeLumpy.py $outDir
chmod u+x $outDir/lumpy/runLumpy.sh
$outDir/lumpy/runLumpy.sh

##make pindel configs
python3 makePindelCfg.py $outDir
now=$(date)
chmod u+x $outDir/runPindel.sh -i $outDir -r $reference -d $now
##TODO nohup $outDir/runPindel.sh -i $outDir -r $reference -d $now


##finish running delly 

delly merge -o $outDir/delly/all_sites $outDir/delly/*.bcf

for b in $bams/*.bam
do
tmp=${b##*/}
dellyOut2=${tmp%.sorted.bam}_delly.bcf
delly call -g $reference -v $outDir/delly/all_site -o $outDir/delly/$dellyOut2 $b
done

