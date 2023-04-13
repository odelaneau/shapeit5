#!/bin/bash

for CHR in 20; do
	ANN=/mnt/project/Phasing/PhasingSNParray/step2_chrrename/chr_rename.txt

	VCFF=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/full_c$CHR\_b0_v2.vcf.gz
	
	OUTF=full_c$CHR\_b0_v2.b37.vcf.gz
	
	dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step2_chrrename/" -icmd="bcftools annotate -Oz -o $OUTF --rename-chrs $ANN $VCFF && bcftools index $OUTF" --tag updatechr --tag chr$CHR --instance-type mem2_ssd1_v2_x2 --name updatechr3_chr$CHR --priority normal -y
	
done
