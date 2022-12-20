#!/bin/bash

## input vcf, output s coefficients and genotype matrix

VCF=$1
SOUTPUT=$2
MOUTPUT=$3
DOUTPUT=$4

module load vital-it/7 UHTS/Analysis/samtools/1.10

bcftools query -f '%INFO/S\n' $VCF > $SOUTPUT

bcftools query -f '%INFO/DOM\n' $VCF > $DOUTPUT

grep -v "^#" $VCF | \
    sed -e 's;\.;NA;g' | \
    awk -F$'\t' '{print substr($0, index($0,$9))}' | \
    awk -F$'\t' '{print substr($0, index($0,$2))}' > $MOUTPUT

# sed -i -e 's;0|0;0;g'  $MOUTPUT
# sed -i -e 's;0|1;0.3;g'  $MOUTPUT # h=0.3
# sed -i -e 's;1|0;0.3;g'  $MOUTPUT # h=0.3
# sed -i -e 's;1|1;2;g'  $MOUTPUT

