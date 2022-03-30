#!/bin/bash

#GET VALIDATION DATA
while read LINE; do
	
	REG=$(echo $LINE | awk '{ print $3; }')
	
	BGL=beagle.25Mar22.4f6.jar
	VCF=/mnt/project/Phasing/PhasingWGS/step2_splitchunks/benchmark_ukb23352_c20_qc_v1.$REG\.vcf.gz
	MAP=/mnt/project/data/shapeit_maps/chr20.b38.gmap.gz
	OUT=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit.default.vcf.gz
	LOG=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit.default.log
	
	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="SHAPEIT5_phase1_static --input $VCF --map $MAP --output $OUT --thread 32 --log $LOG --filter-maf 0.001 --filter-snp && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_shapeit5_$REG -y
	
done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt
