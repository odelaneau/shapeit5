#!/bin/bash

VCF=/mnt/project/Phasing/PhasingSNParray/step4_liftover/full_c20_b0_v2.b38.sorted.vcf.gz

#GET VALIDATION
SUB=/mnt/project/Phasing/PhasingSNParray/step5_benchmark/validation.samples.txt
OUT=validation_c20_b0_v2.b38.sorted.vcf.gz
dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step5_benchmark/" -icmd="bcftools view -Oz -o $OUT -S $SUB --force-samples $VCF && bcftools index $OUT" --tag subset --tag benchmark --instance-type mem2_ssd1_v2_x2 --name benchmark_subset --priority normal -y

#GET DOWNSAMPLES
for N in 5000 10000 20000 50000 100000 200000 300000 400000 480853; do
	OUT=benchmark_c20_b0_v2.b38.sorted.N$N\.vcf.gz
	SUB=/mnt/project/Phasing/PhasingSNParray/step5_benchmark/N$N\.txt
	dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step5_benchmark/" -icmd="bcftools view -Oz -o $OUT -S $SUB $VCF && bcftools index $OUT" --tag subset --tag benchmark --instance-type mem2_ssd1_v2_x2 --name benchmark_subset --priority normal -y
done

#PHASING
for N in 5000 10000 20000 50000 100000 200000 300000 400000 480853; do
	MAP=/mnt/project/data/shapeit_maps/chr20.b38.gmap.gz
	INP=/mnt/project/Phasing/PhasingSNParray/step5_benchmark/benchmark_c20_b0_v2.b38.sorted.N$N\.bcf
		
	LOG=benchmark_c20_b0_v2.b38.sorted.N$N\.phased.log
	OUT=benchmark_c20_b0_v2.b38.sorted.N$N\.phased.bcf

	dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step5_benchmark/" -iimage_file="/docker/shapeit5_0.0.1.tar.gz" -icmd="SHAPEIT5_phase1_static --input $INP --map $MAP --output $OUT --region chr$CHR --log $LOG --thread 16 && bcftools index $OUT --threads 16" --tag phasing --tag benchmark --instance-type mem2_ssd1_v2_x16 --priority normal --name benchmark_phasing -y
done

