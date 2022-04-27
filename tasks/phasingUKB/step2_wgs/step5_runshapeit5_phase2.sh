#!/bin/bash

while read LINE; do
	SREG=$(echo $LINE | awk '{ print $3; }')
	IREG=$(echo $LINE | awk '{ print $4; }')
	BCF=/mnt/project/Phasing/PhasingWGS/step2_splitchunks/benchmark_ukb23352_c20_qc_v1.$SREG\.bcf
	SCA1=/mnt/project/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.indel.default.depth8.bcf
	SCA2=/mnt/project/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.indel.scaffold.depth8.bcf
	MAP=/mnt/project/data/shapeit_maps/chr20.b38.gmap.gz
	
	OUT=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.default.bcf
	LOG=$(basename $OUT .bcf)\.log
	TIM=$(basename $OUT .bcf)\.time
	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step6_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase2_static --input $BCF --scaffold $SCA1 --map $MAP --output $OUT --log $LOG --scaffold-region $SREG --input-region $IREG --thread 32 && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag phase2 --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_shapeit5a_$SREG -y

	OUT=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.scaffold.bcf
	LOG=$(basename $OUT .bcf)\.log
	TIM=$(basename $OUT .bcf)\.time
	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step6_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase2_static --input $BCF --scaffold $SCA2 --map $MAP --output $OUT --log $LOG --scaffold-region $SREG --input-region $IREG --thread 32 && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag phase2 --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_shapeit5b_$SREG -y

done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt

