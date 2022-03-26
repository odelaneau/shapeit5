#!/bin/bash

for CHR in 20; do
	ANN=/mnt/project/Phasing/PhasingSNParray/step2_chrrename/chr_rename.txt

	VCFB=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/benchmark_c$CHR\_b0_v2.vcf.gz
	VCFF=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/full_c$CHR\_b0_v2.vcf.gz
	VCFV=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/validation_c$CHR\_b0_v2.vcf.gz
	
	OUTB=benchmark_c$CHR\_b0_v2.b37.vcf.gz
	OUTF=full_c$CHR\_b0_v2.b37.vcf.gz
	OUTV=validation_c$CHR\_b0_v2.b37.vcf.gz
	
	dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step2_chrrename/" -icmd="bcftools annotate -Oz -o $OUTB --rename-chrs $ANN $VCFB && bcftools index $OUTB" --tag updatechr --tag chr$CHR --instance-type mem2_ssd1_v2_x2 --name updatechr1_chr$CHR --priority normal -y
	
	dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step2_chrrename/" -icmd="bcftools annotate -Oz -o $OUTV --rename-chrs $ANN $VCFV && bcftools index $OUTV" --tag updatechr --tag chr$CHR --instance-type mem2_ssd1_v2_x2 --name updatechr2_chr$CHR --priority normal -y
		
	dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step2_chrrename/" -icmd="bcftools annotate -Oz -o $OUTF --rename-chrs $ANN $VCFF && bcftools index $OUTF" --tag updatechr --tag chr$CHR --instance-type mem2_ssd1_v2_x2 --name updatechr3_chr$CHR --priority normal -y
	
done
