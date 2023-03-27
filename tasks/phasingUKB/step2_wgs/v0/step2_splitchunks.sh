#!/bin/bash

for N in 2000 5000 10000 20000 50000 100000 147754; do
	BCF=/mnt/project/Phasing/PhasingWGS/step1_preparedata/benchmark_ukb23352_c20_qc_v1.subset.N$N\.bcf
	
	while read LINE; do
		REG=$(echo $LINE | awk '{ print $3; }')
		OUTV=$(basename $BCF .bcf)\.$REG\.vcf.gz
		OUTB=$(basename $BCF .bcf)\.$REG\.bcf
		
		dx run app-swiss-army-knife --folder "/Phasing/PhasingWGS/step2_splitchunks/" -icmd="bcftools view -Oz -o $OUTV $BCF $REG && bcftools index $OUTV" --tag splitB --tag benchWGS --instance-type mem2_ssd1_v2_x2 --name benchWGS_splitB --priority normal -y
		dx run app-swiss-army-knife --folder "/Phasing/PhasingWGS/step2_splitchunks/" -icmd="bcftools view -Ob -o $OUTB $BCF $REG && bcftools index $OUTB" --tag splitB --tag benchWGS --instance-type mem2_ssd1_v2_x2 --name benchWGS_splitB --priority normal -y
		
	done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt
done

