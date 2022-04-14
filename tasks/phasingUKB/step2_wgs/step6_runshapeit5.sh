#!/bin/bash

# WITHOUT PIRS
while read LINE; do

	SREG=$(echo $LINE | awk '{ print $3; }')
	IREG=$(echo $LINE | awk '{ print $4; }')
	BCF=/mnt/project/Phasing/PhasingWGS/step2_splitchunks/benchmark_ukb23352_c20_qc_v1.$SREG\.bcf
	SCA=/mnt/project/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.default.depth8.bcf
	MAP=/mnt/project/data/shapeit_maps/chr20.b38.gmap.gz
	
	OUT=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.phase2.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.phase2.log
	TIM=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.phase2.time
	
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step6_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase2_static --input $BCF --scaffold $SCA --map $MAP --output $OUT --log $LOG --scaffold-region $SREG --input-region $IREG --thread 32 && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5.2 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority low --name benchWGS_shapeit5.2_$SREG -y
done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt

