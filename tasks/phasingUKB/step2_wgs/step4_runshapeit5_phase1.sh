#!/bin/bash

MAP=/mnt/project/data/shapeit_maps/chr20.b38.gmap.gz
SCA=/mnt/project/Phasing/PhasingSNParray/step5_benchmark/benchmark_c20_b0_v2.b38.sorted.N480853.phased.bcf

DOCKER=shapeit5_$(git log -1 --format=%cd --date=short)\_$(git rev-parse --short HEAD)\.tar.gz

for N in 2000 5000 10000 20000 50000 100000 147754; do

	while read LINE; do
		
		REG=$(echo $LINE | awk '{ print $3; }')
		BCF=/mnt/project/Phasing/PhasingWGS/step2_splitchunks/chunks/benchmark_ukb23352_c20_qc_v1.subset.N$N\.$REG\.bcf
		
		OUT=benchmark_ukb23352_c20_qc_v1.subset.N$N\.$REG\.shapeit5.default.bcf
		LOG=benchmark_ukb23352_c20_qc_v1.subset.N$N\.$REG\.shapeit5.default.log
		TIM=benchmark_ukb23352_c20_qc_v1.subset.N$N\.$REG\.shapeit5.default.time
		
		#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step4_runshapeit_phase1/N$N/chunks/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase_common_static --input $BCF --map $MAP --output $OUT --thread 32 --log $LOG --filter-maf 0.001 --region $REG && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority high --name benchWGS_shapeit5_$REG -y
		
	done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt
done




for N in 2000 5000 10000 20000 50000 100000 147754; do
	#Produce list of files
	rm benchmark_ukb23352_c20_qc_v1.subset.N$N\.shapeit5.default.list
	while read LINE; do
		REG=$(echo $LINE | awk '{ print $3; }')
		echo /mnt/project/Phasing/PhasingWGS/step4_runshapeit_phase1/N$N/chunks/benchmark_ukb23352_c20_qc_v1.subset.N$N\.$REG\.shapeit5.default.bcf >> benchmark_ukb23352_c20_qc_v1.subset.N$N\.shapeit5.default.list
	done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt
	dx cd /Phasing/PhasingWGS/step4_runshapeit_phase1/N$N
	dx rm benchmark_ukb23352_c20_qc_v1.subset.N$N\.shapeit5.default.list
	dx upload benchmark_ukb23352_c20_qc_v1.subset.N$N\.shapeit5.default.list
	rm benchmark_ukb23352_c20_qc_v1.subset.N$N\.shapeit5.default.list
	
	#RUN ligate
	TXT=/mnt/project/Phasing/PhasingWGS/step4_runshapeit_phase1/N$N/benchmark_ukb23352_c20_qc_v1.subset.N$N\.shapeit5.default.list
	TIM=benchmark_ukb23352_c20_qc_v1.subset.N$N\.fullchr.shapeit5.default.time
	OUT=benchmark_ukb23352_c20_qc_v1.subset.N$N\.fullchr.shapeit5.default.bcf
		
	dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step4_runshapeit_phase1/N$N/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_ligate_static --input $TXT --output $OUT --thread 2 --index" --tag benchWGS --tag ligate --tag chr20 --instance-type mem2_ssd1_v2_x2 --priority normal --name benchWGS_ligate -y
done
	
	 
	
	





