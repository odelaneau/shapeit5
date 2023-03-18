#!/bin/bash

MAP=/mnt/project/data/shapeit_maps/chr20.b38.gmap.gz

DOCKER=shapeit5_$(git log -1 --format=%cd --date=short)\_$(git rev-parse --short HEAD)\.tar.gz

for N in 2000 5000 10000 20000 50000 100000 147754; do
	while read LINE; do
		
		SREG=$(echo $LINE | awk '{ print $3; }')
		IREG=$(echo $LINE | awk '{ print $4; }')
		BCF=/mnt/project/Phasing/PhasingWGS/step1_preparedata/benchmark_ukb23352_c20_qc_v1.subset.N$N\.bcf
		
		SCA=/mnt/project/Phasing/PhasingWGS/step4_runshapeit_phase1/N$N\.test/benchmark_ukb23352_c20_qc_v1.subset.N$N\.fullchr.shapeit5.default.bcf
		OUT=benchmark_ukb23352_c20_qc_v1.subset.N$N\.$SREG\.shapeit5.default.bcf
		TIM=benchmark_ukb23352_c20_qc_v1.subset.N$N\.$SREG\.shapeit5.default.time
		LOG=benchmark_ukb23352_c20_qc_v1.subset.N$N\.$SREG\.shapeit5.default.log
				
		#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase_rare_static --input-plain $BCF --scaffold $SCA --map $MAP --output $OUT --log $LOG --scaffold-region $SREG --input-region $IREG --thread 32 && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag phase2 --tag $SREG --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_shapeit5_phase2a_$SREG -y

	done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/v0/step2_splitchunks/chr20.size4Mb.txt

done




DOCKER=shapeit5.test.tar.gz

while read LINE; do
	SREG=$(echo $LINE | awk '{ print $3; }')
	IREG=$(echo $LINE | awk '{ print $4; }')
	
	#PCF=/mnt/project/data/ukb_wgs/unphased/qc/ukb23352_c20_qc_v1.bcf
	#PED=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/samples.families.caucasian.txt
	#SCA=/mnt/project/Phasing/PhasingWGS/step9_production/ukb23352_c20_qc_v1.common.phased.bcf
	#SCA=/mnt/project/Phasing/PhasingWGS/step4_runshapeit_phase1/N$N\.test/benchmark_ukb23352_c20_qc_v1.subset.N147754.fullchr.shapeit5.ped.bcf
	BCF=/mnt/project/Phasing/PhasingWGS/step1_preparedata/benchmark_ukb23352_c20_qc_v1.subset.N147754.bcf
	SCA=/mnt/project/Phasing/PhasingWGS/step4_runshapeit_phase1/N$N\.test/benchmark_ukb23352_c20_qc_v1.subset.N147754.fullchr.shapeit5.default.bcf
	
	OUT=benchmark_ukb23352_c20_qc_v1.subset.N147754.$SREG\.shapeit5.default.bcf
	TIM=benchmark_ukb23352_c20_qc_v1.subset.N147754.$SREG\.shapeit5.default.time
	LOG=benchmark_ukb23352_c20_qc_v1.subset.N147754.$SREG\.shapeit5.default.log
	
	dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N147754.test/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase_rare_static --input $BCF --scaffold $SCA --map $MAP --output $OUT --log $LOG --scaffold-region $SREG --input-region $IREG --thread 16" --tag benchWGS --tag shapeit5 --tag phase2 --tag $SREG --instance-type mem3_ssd1_v2_x16 --priority normal --name benchWGS_shapeit5_phase2a_$SREG -y

done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/v0/step2_splitchunks/chr20.size4Mb.txt

