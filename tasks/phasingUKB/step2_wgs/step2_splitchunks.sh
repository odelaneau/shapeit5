#!/bin/bash

BCF=/mnt/project/Phasing/PhasingWGS/step1_preparedata/benchmark_ukb23352_c20_qc_v1.bcf

#GET VALIDATION DATA
while read LINE; do
	
	REG=$(echo $LINE | awk '{ print $3; }')
	OUT=$(basename $BCF .bcf)\.$REG\.vcf.gz
	
	dx run app-swiss-army-knife --folder "/Phasing/PhasingWGS/step2_splitchunks/" -icmd="bcftools view -Oz -o $OUT $BCF $REG && bcftools index $OUT" --tag splitB --tag benchWGS --instance-type mem2_ssd1_v2_x2 --name benchWGS_splitB --priority normal -y
done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt
