#!/bin/bash

for CHR in 20; do
	MAP=/mnt/project/data/shapeit_maps/chr$CHR\.b38.gmap.gz
	VCF=/mnt/project/Phasing/PhasingSNParray/step4_liftover/benchmark_c$CHR\_b0_v2.b38.sorted.vcf.gz
	
	LOG=benchmark_c$CHR\_b0_v2.b38.sorted.log
	BCF=benchmark_c$CHR\_b0_v2.b38.sorted.bcf

	dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step5_phasing/" -iimage_file="/docker/shapeit5_0.0.1.tar.gz" -icmd="SHAPEIT5_phase1_static --input $VCF --map $MAP --output $BCF --region chr$CHR --log $LOG --thread 32 && bcftools index $BCF --threads 32" --tag phasingSA --tag chr20 --instance-type mem2_ssd1_v2_x32 --priority normal --name phasingSA_chr$CHR -y

done

