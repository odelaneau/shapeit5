#!/bin/bash

MAP=/mnt/project/data/shapeit_maps/chr20.b38.gmap.gz
SCA=/mnt/project/Phasing/PhasingSNParray/step5_benchmark/benchmark_c20_b0_v2.b38.sorted.N480853.phased.bcf
PED=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/samples.families.caucasian.txt

#DOCKER=shapeit5_$(git log -1 --format=%cd --date=short)\_$(git rev-parse --short HEAD)\.tar.gz
DOCKER=shapeit5.test.tar.gz

#for N in 2000 5000 10000 20000 50000 100000 147754; do
for N in 147754; do

	#Produce list of files
	rm benchmark_ukb23352_c20_qc_v1.subset.N147754\.shapeit5.ped.list
	while read LINE; do
		REG=$(echo $LINE | awk '{ print $3; }')
		echo /mnt/project/Phasing/PhasingWGS/step4_runshapeit_phase1/N147754.test/chunks/benchmark_ukb23352_c20_qc_v1.subset.N147754.$REG\.shapeit5.ped.bcf >> benchmark_ukb23352_c20_qc_v1.subset.N147754\.shapeit5.ped.list
	done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/v0/step2_splitchunks/chr20.size4Mb.txt
	dx cd /Phasing/PhasingWGS/step4_runshapeit_phase1/N147754.test
	dx rm benchmark_ukb23352_c20_qc_v1.subset.N147754.shapeit5.ped.list
	dx upload benchmark_ukb23352_c20_qc_v1.subset.N147754.shapeit5.ped.list
	rm benchmark_ukb23352_c20_qc_v1.subset.N147754.shapeit5.ped.list
	
	#RUN ligate
	TXT=/mnt/project/Phasing/PhasingWGS/step4_runshapeit_phase1/N147754.test/benchmark_ukb23352_c20_qc_v1.subset.N147754.shapeit5.ped.list
	TIM=benchmark_ukb23352_c20_qc_v1.subset.N147754.fullchr.shapeit5.ped.time
	OUT=benchmark_ukb23352_c20_qc_v1.subset.N147754.fullchr.shapeit5.ped.bcf
		
	dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step4_runshapeit_phase1/N147754.test/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_ligate_static --input $TXT --output $OUT --pedigree $PED --thread 2 --index" --tag benchWGS --tag ligate --tag chr20 --instance-type mem2_ssd1_v2_x2 --priority normal --name benchWGS_ligate -y
done
