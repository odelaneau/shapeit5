#!/bin/bash

MAP=/mnt/project/data/shapeit_maps/chr20.b38.gmap.gz

DOCKER=shapeit5_$(git log -1 --format=%cd --date=short)\_$(git rev-parse --short HEAD)\.tar.gz



for N in 2000 5000 10000 20000 50000 100000 147754; do

	while read LINE; do
		
		SREG=$(echo $LINE | awk '{ print $3; }')
		IREG=$(echo $LINE | awk '{ print $4; }')
		BCF=/mnt/project/Phasing/PhasingWGS/step2_splitchunks/benchmark_ukb23352_c20_qc_v1.subset.N$N\.$SREG\.bcf

		#NOT SCAFFOLDED		
		SCA=/mnt/project/Phasing/PhasingWGS/step4_runshapeit_phase1/N$N/benchmark_ukb23352_c20_qc_v1.subset.N$N\.$SREG\.shapeit5.default.bcf
		OUT=benchmark_ukb23352_c20_qc_v1.subset.N$N\.$SREG\.shapeit5.default.bcf
		TIM=benchmark_ukb23352_c20_qc_v1.subset.N$N\.$SREG\.shapeit5.default.time
		LOG=benchmark_ukb23352_c20_qc_v1.subset.N$N\.$SREG\.shapeit5.default.log
				
		#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase_rare_static --input-plain $BCF --scaffold $SCA --map $MAP --output $OUT --log $LOG --scaffold-region $SREG --input-region $IREG --thread 32 && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag phase2 --tag $SREG --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_shapeit5_phase2a_$SREG -y
		
	done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt

done


for N in 5000 10000 20000 50000 100000 147754; do

	while read LINE; do
		
		SREG=$(echo $LINE | awk '{ print $3; }')
		IREG=$(echo $LINE | awk '{ print $4; }')
		BCF=/mnt/project/Phasing/PhasingWGS/step1_preparedata/benchmark_ukb23352_c20_qc_v1.subset.N$N\.bcf
		
		SCA=/mnt/project/Phasing/PhasingWGS/step4_runshapeit_phase1/N$N/benchmark_ukb23352_c20_qc_v1.subset.N$N\.fullchr.shapeit5.default.bcf
		OUT=benchmark_ukb23352_c20_qc_v1.subset.N$N\.$SREG\.shapeit5.ligated.bcf
		TIM=benchmark_ukb23352_c20_qc_v1.subset.N$N\.$SREG\.shapeit5.ligated.time
		LOG=benchmark_ukb23352_c20_qc_v1.subset.N$N\.$SREG\.shapeit5.ligated.log
				
		dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase_rare_static --input-plain $BCF --scaffold $SCA --map $MAP --output $OUT --log $LOG --scaffold-region $SREG --input-region $IREG --thread 32 && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag phase2 --tag $SREG --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_shapeit5_phase2a_$SREG -y
		
	done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt

done




