#!/bin/bash

#GET VALIDATION DATA
while read LINE; do

	SREG=$(echo $LINE | awk '{ print $3; }')
	IREG=$(echo $LINE | awk '{ print $4; }')
	BCF=/mnt/project/Phasing/PhasingWGS/step2_splitchunks/benchmark_ukb23352_c20_qc_v1.$IREG\.bcf
	SCA=/mnt/project/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.$IREG\.shapeit5.default.depth8.bcf
	MAP=/mnt/project/data/shapeit_maps/chr20.b38.gmap.gz
	
	OUT=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.default.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.default.log
	TIM=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.default.time
	
	
	

	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase1_static --input $BCF --map $MAP --output $OUT --thread 32 --log $LOG --filter-maf 0.001 --filter-snp --region $REG && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority low --name benchWGS_shapeit5_$REG -y

	OUT=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.scaffold.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.scaffold.log
	TIM=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.scaffold.time

#	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase1_static --input $BCF --scaffold $SCA --map $MAP --output $OUT --thread 32 --log $LOG --filter-maf 0.001 --filter-snp --region $REG && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority low --name benchWGS_shapeit5_$REG -y
	
	OUT=
	LOG=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.default.depth8.log
	TIM=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.default.depth8.time

	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase1_static --input $BCF --map $MAP --output $OUT --thread 32 --log $LOG --filter-maf 0.001 --filter-snp --pbwt-depth 8 --region $REG && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority low --name benchWGS_shapeit5_$REG -y

	OUT=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.scaffold.depth8.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.scaffold.depth8.log
	TIM=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.scaffold.depth8.time

	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase1_static --input $BCF --scaffold $SCA --map $MAP --output $OUT --thread 32 --log $LOG --filter-maf 0.001 --filter-snp --pbwt-depth 8 --region $REG && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority low --name benchWGS_shapeit5_$REG -y

	OUT=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.default.modulo5.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.default.modulo5.log
	TIM=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.default.modulo5.time

	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase1_static --input $BCF --map $MAP --output $OUT --thread 32 --log $LOG --filter-maf 0.001 --filter-snp --pbwt-modulo 0.05 --region $REG && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority low --name benchWGS_shapeit5_$REG -y

	OUT=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.scaffold.modulo5.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.scaffold.modulo5.log
	TIM=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.scaffold.modulo5.time

	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase1_static --input $BCF --scaffold $SCA --map $MAP --output $OUT --thread 32 --log $LOG --filter-maf 0.001 --filter-snp --pbwt-modulo 0.05 --region $REG && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority low --name benchWGS_shapeit5_$REG -y
	
	OUT=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.scaffold.depth16.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.scaffold.depth16.log
	TIM=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.scaffold.depth16.time

	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase1_static --input $BCF --scaffold $SCA --map $MAP --output $OUT --thread 32 --log $LOG --filter-maf 0.001 --filter-snp --pbwt-depth 16 --region $REG && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority low --name benchWGS_shapeit5_$REG -y


done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt


