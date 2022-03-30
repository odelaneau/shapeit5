#!/bin/bash

#GET VALIDATION DATA
while read LINE; do
	
	REG=$(echo $LINE | awk '{ print $3; }')
	
	BGL=beagle.25Mar22.4f6.jar
	VCF=/mnt/project/Phasing/PhasingWGS/step2_splitchunks/benchmark_ukb23352_c20_qc_v1.$REG\.vcf.gz
	MAP=/mnt/project/data/plink_maps/plink.prefix.chr20.GRCh38.map
	OUT=benchmark_ukb23352_c20_qc_v1.$REG\.beagle5.3.vcf.gz
	TIM=benchmark_ukb23352_c20_qc_v1.$REG\.beagle5.3.time
	
	dx run app-swiss-army-knife -iin="/docker/beagle.25Mar22.4f6.jar" --folder="/Phasing/PhasingWGS/step3_runbeagle/" -icmd="/usr/bin/time -vo $TIM java -Xmx256G -jar $BGL gt=$VCF map=$MAP out=$OUT window=500.0 nthreads=32 chrom=chr20 && bcftools index -f $OUT --threads 8" --tag benchWGS --tag beagle5.3 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_beagle5.3_$REG -y
done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt




		
