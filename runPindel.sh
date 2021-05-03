#!/bin/bash

while getopts i:r:d flag
do
    case "${flag}" in
        i) pinDir=${OPTARG};;
        r) reference=${OPTARG};;
        d) date=${OPTARG};;
    esac
done

cd $outDir/pindel

for f in *.cfg
do
tmp=${f##*/}
sammy=${tmp%_pindel.cfg}
pVcf=${f%.cfg}.vcf
pindel -f $reference -i $f -o $sammy -c ALL -T 8 -a 4 -d 40 -m 4 
/usr/local/pindel/pindel2vcf -P $sammy -r $reference -R $reference -d $date -v $pVcf
done



