#!/bin/bash

while read LINE; do
	SREG=$(echo $LINE | awk '{ print $3; }')
	IREG=$(echo $LINE | awk '{ print $4; }')
	BCF=/mnt/project/Phasing/PhasingWGS/step2_splitchunks/benchmark_ukb23352_c20_qc_v1.$SREG\.bcf
	SCA=/mnt/project/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.default.depth8.bcf
	MAP=/mnt/project/data/shapeit_maps/chr20.b38.gmap.gz

	OUT=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.dc8.dr2.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.dc8.dr2.log
	TIM=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.dc8.dr2.time
	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step6_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase2_static --input $BCF --scaffold $SCA --map $MAP --output $OUT --log $LOG --scaffold-region $SREG --input-region $IREG --thread 32 --pbwt-depth-common 8 --pbwt-depth-rare 2 && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5.2 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_shapeit5a_$SREG -y

	OUT=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.dc4.dr2.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.dc4.dr2.log
	TIM=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.dc4.dr2.time
	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step6_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase2_static --input $BCF --scaffold $SCA --map $MAP --output $OUT --log $LOG --scaffold-region $SREG --input-region $IREG --thread 32 --pbwt-depth-common 4 --pbwt-depth-rare 2 && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5.2 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_shapeit5b_$SREG -y

	OUT=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.dc2.dr2.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.dc2.dr2.log
	TIM=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.dc2.dr2.time
	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step6_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase2_static --input $BCF --scaffold $SCA --map $MAP --output $OUT --log $LOG --scaffold-region $SREG --input-region $IREG --thread 32 --pbwt-depth-common 2 --pbwt-depth-rare 2 && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5.2 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_shapeit5c_$SREG -y

done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt



while read LINE; do
	SREG=$(echo $LINE | awk '{ print $3; }')
	IREG=$(echo $LINE | awk '{ print $4; }')
	BCF=/mnt/project/Phasing/PhasingWGS/step2_splitchunks/benchmark_ukb23352_c20_qc_v1.$SREG\.bcf
	SCA=/mnt/project/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.default.depth8.bcf
	MAP=/mnt/project/data/shapeit_maps/chr20.b38.gmap.gz

	OUT=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.b2.m5.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.b2.m5.log
	TIM=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.b2.m5.time
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step6_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase2_static --input $BCF --scaffold $SCA --map $MAP --output $OUT --log $LOG --scaffold-region $SREG --input-region $IREG --thread 32 --mcmc-iterations 5 --mcmc-burnin 2 && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5.2 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority low --name benchWGS_shapeit5a_$SREG -y

	OUT=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.b5.m10.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.b5.m10.log
	TIM=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.b5.m10.time
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step6_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase2_static --input $BCF --scaffold $SCA --map $MAP --output $OUT --log $LOG --scaffold-region $SREG --input-region $IREG --thread 32 --mcmc-iterations 10 --mcmc-burnin 5 && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5.2 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority low --name benchWGS_shapeit5b_$SREG -y

	OUT=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.b5.m20.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.b5.m20.log
	TIM=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.b5.m20.time
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step6_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase2_static --input $BCF --scaffold $SCA --map $MAP --output $OUT --log $LOG --scaffold-region $SREG --input-region $IREG --thread 32 --mcmc-iterations 20 --mcmc-burnin 5 && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5.2 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority low --name benchWGS_shapeit5c_$SREG -y

	OUT=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.b10.m20.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.b10.m20.log
	TIM=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.b10.m20.time
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step6_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase2_static --input $BCF --scaffold $SCA --map $MAP --output $OUT --log $LOG --scaffold-region $SREG --input-region $IREG --thread 32 --mcmc-iterations 20 --mcmc-burnin 10 && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5.2 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority low --name benchWGS_shapeit5d_$SREG -y

	OUT=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.b10.m50.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.b10.m50.log
	TIM=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.b10.m50.time
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step6_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase2_static --input $BCF --scaffold $SCA --map $MAP --output $OUT --log $LOG --scaffold-region $SREG --input-region $IREG --thread 32 --mcmc-iterations 50 --mcmc-burnin 10 && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5.2 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority low --name benchWGS_shapeit5e_$SREG -y

done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt

