#!/bin/bash

BCF=/mnt/project/data/ukb_wgs/unphased/qc/ukb23352_c20_qc_v1.bcf

#GET VALIDATION DATA
SUB=/mnt/project/Phasing/PhasingWGS/step1_preparedata/validation.samples.txt
OUT=validation_ukb23352_c20_qc_v1.bcf
#dx run app-swiss-army-knife --folder "/Phasing/PhasingWGS/step1_preparedata/" -icmd="bcftools view -Ob -o $OUT -S $SUB --force-samples $BCF && bcftools index $OUT" --tag subsetV --tag benchWGS --instance-type mem2_ssd1_v2_x2 --name benchWGS_subsetV --priority normal -y

#GET BENCHMARK DATA
for N in 2000 5000 10000 20000 50000 100000 147754; do
	SUB=/mnt/project/Phasing/PhasingWGS/step1_preparedata/N$N\.merged.txt
	OUT=benchmark_ukb23352_c20_qc_v1.subset.N$N\.bcf
	#dx run app-swiss-army-knife --folder "/Phasing/PhasingWGS/step1_preparedata/" -icmd="bcftools view -Ou -S $SUB --force-samples $BCF | bcftools view -Ob -o $OUT -c 1:minor - && bcftools index $OUT" --tag subsetN --tag benchWGS --instance-type mem2_ssd1_v2_x2 --name benchWGS_sub$N --priority normal -y
done

#GET SITE LISTS
for N in 2000 5000 10000 20000 50000 100000 147754; do
	BCF=/mnt/project/Phasing/PhasingWGS/step1_preparedata/benchmark_ukb23352_c20_qc_v1.subset.N$N\.bcf
	FQC=frequencies.subset.N$N\.SNPs.common.bcf
	FQR=frequencies.subset.N$N\.SNPs.all.bcf
	FQF=frequencies.subset.N$N\.ALL.all.bcf

	dx run app-swiss-army-knife --folder "/Phasing/PhasingWGS/step1_preparedata/" -icmd="bcftools view -G --min-af 0.001:minor -m2 -M2 -v snps -Ob -o $FQC $BCF && bcftools index $FQC" --tag subsetN --tag benchWGS --instance-type mem2_ssd1_v2_x2 --name benchWGS_sub1 --priority normal -y
	dx run app-swiss-army-knife --folder "/Phasing/PhasingWGS/step1_preparedata/" -icmd="bcftools view -G -m2 -M2 -v snps -Ob -o $FQR $BCF && bcftools index $FQR" --tag subsetN --tag benchWGS --instance-type mem2_ssd1_v2_x2 --name benchWGS_sub2 --priority normal -y
	dx run app-swiss-army-knife --folder "/Phasing/PhasingWGS/step1_preparedata/" -icmd="bcftools view -G -Ob -o $FQF $BCF && bcftools index $FQF" --tag subsetN --tag benchWGS --instance-type mem2_ssd1_v2_x2 --name benchWGS_sub3 --priority normal -y
	
done


#Whole CHR convertion BCF to VCF for Beagle
for N in 2000 5000 10000 20000 50000 100000; do
	BCF=/mnt/project/Phasing/PhasingWGS/step1_preparedata/benchmark_ukb23352_c20_qc_v1.subset.N$N\.bcf
	VCF=$(basename $BCF .bcf)\.vcf.gz
	dx run app-swiss-army-knife --folder "/Phasing/PhasingWGS/step2_splitchunks/" -icmd="bcftools view -Oz -o $VCF $BCF --threads 8 && bcftools index $VCF --threads 8" --tag bcf2vcf --tag benchWGS --instance-type mem2_ssd1_v2_x8 --name benchWGS_conv --priority normal -y
done

