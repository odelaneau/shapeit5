#!/bin/bash

for CHR in 20; do
	MAP=/mnt/project/data/shapeit_maps/chr$CHR\.b37.gmap.gz
	VCF=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/benchmark_c$CHR\_b0_v2.vcf.gz
	
	LOG=benchmark_c$CHR\_b0_v2.phased.log
	BCF=benchmark_c$CHR\_b0_v2.phased.bcf

	dx run app-swiss-army-knife	--folder "/mnt/project/Phasing/PhasingSNParray/step2_phasing/" -iimage_file="/docker/shapeit5\:0.0.1.tar.gz" -icmd="SHAPEIT5_phase1_static --input $VCF --map $MAP --output $BCF --region $CHR --log $LOG --thread 32 && bcftools index -f $OUT --threads 32" --tag phasingSA --tag chr20 --instance-type mem2_ssd1_v2_x32 --priority normal --name shapeit5_phase -y
done