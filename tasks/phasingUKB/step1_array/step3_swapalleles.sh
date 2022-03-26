#!/bin/bash

for CHR in 20; do
	VCFB=/mnt/project/Phasing/PhasingSNParray/step2_chrrename/benchmark_c$CHR\_b0_v2.b37.vcf.gz
	VCFF=/mnt/project/Phasing/PhasingSNParray/step2_chrrename/full_c$CHR\_b0_v2.b37.vcf.gz
	VCFV=/mnt/project/Phasing/PhasingSNParray/step2_chrrename/validation_c$CHR\_b0_v2.b37.vcf.gz

	OUTB=benchmark_c$CHR\_b0_v2.b37.swapped.vcf.gz
	OUTF=full_c$CHR\_b0_v2.b37.swapped.vcf.gz
	OUTV=validation_c$CHR\_b0_v2.b37.swapped.vcf.gz
	
	dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step3_swapalleles/" -iimage_file="/Phasing/PhasingSNParray/step3_swapalleles/swaprefalt_0.0.1.tar.gz" -icmd="swapRefAlt_static --input $VCFV --output $OUTV && bcftools index $OUTV" --tag swapalleles --tag chr$CHR --instance-type mem2_ssd1_v2_x2 --priority normal --name swapalleles1_chr$CHR -y
	
	dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step3_swapalleles/" -iimage_file="/Phasing/PhasingSNParray/step3_swapalleles/swaprefalt_0.0.1.tar.gz" -icmd="swapRefAlt_static --input $VCFF --output $OUTF && bcftools index $OUTF" --tag swapalleles --tag chr$CHR --instance-type mem2_ssd1_v2_x2 --priority normal --name swapalleles2_chr$CHR -y
	
	dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step3_swapalleles/" -iimage_file="/Phasing/PhasingSNParray/step3_swapalleles/swaprefalt_0.0.1.tar.gz" -icmd="swapRefAlt_static --input $VCFB --output $OUTB && bcftools index $OUTB" --tag swapalleles --tag chr$CHR --instance-type mem2_ssd1_v2_x2 --priority normal --name swapalleles3_chr$CHR -y
done
