#!/bin/bash

#GET VALIDATION DATA
while read LINE; do

	REG=$(echo $LINE | awk '{ print $3; }')
	BCF=/mnt/project/Phasing/PhasingWGS/step2_splitchunks/benchmark_ukb23352_c20_qc_v1.$REG\.bcf
	SCA=/mnt/project/Phasing/PhasingSNParray/step5_benchmark/benchmark_c20_b0_v2.b38.sorted.N480853.phased.bcf
	MAP=/mnt/project/data/shapeit_maps/chr20.b38.gmap.gz

	OUT=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.default.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.default.log
	TIM=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.default.time

	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase1_static --input $BCF --map $MAP --output $OUT --thread 32 --log $LOG --filter-maf 0.001 --region $REG && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_shapeit5_$REG -y

	OUT=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.scaffold.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.scaffold.log
	TIM=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.scaffold.time

	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase1_static --input $BCF --scaffold $SCA --map $MAP --output $OUT --thread 32 --log $LOG --filter-maf 0.001 --region $REG && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_shapeit5_$REG -y
	
	OUT=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.default.depth8.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.default.depth8.log
	TIM=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.default.depth8.time

	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase1_static --input $BCF --map $MAP --output $OUT --thread 32 --log $LOG --filter-maf 0.001 --pbwt-depth 8 --region $REG && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_shapeit5_$REG -y

	OUT=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.scaffold.depth8.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.scaffold.depth8.log
	TIM=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.scaffold.depth8.time

	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase1_static --input $BCF --scaffold $SCA --map $MAP --output $OUT --thread 32 --log $LOG --filter-maf 0.001 --pbwt-depth 8 --region $REG && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_shapeit5_$REG -y

	OUT=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.default.modulo5.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.default.modulo5.log
	TIM=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.default.modulo5.time

	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase1_static --input $BCF --map $MAP --output $OUT --thread 32 --log $LOG --filter-maf 0.001 --pbwt-modulo 0.05 --region $REG && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_shapeit5_$REG -y

	OUT=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.scaffold.modulo5.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.scaffold.modulo5.log
	TIM=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.scaffold.modulo5.time

	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase1_static --input $BCF --scaffold $SCA --map $MAP --output $OUT --thread 32 --log $LOG --filter-maf 0.001 --pbwt-modulo 0.05 --region $REG && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_shapeit5_$REG -y

	OUT=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.default.modulo25.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.default.modulo25.log
	TIM=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.default.modulo25.time

	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase1_static --input $BCF --map $MAP --output $OUT --thread 32 --log $LOG --filter-maf 0.001 --pbwt-modulo 0.025 --region $REG && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_shapeit5_$REG -y

	OUT=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.scaffold.modulo25.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.scaffold.modulo25.log
	TIM=benchmark_ukb23352_c20_qc_v1.$REG\.shapeit5.indel.scaffold.modulo25.time

	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase1_static --input $BCF --scaffold $SCA --map $MAP --output $OUT --thread 32 --log $LOG --filter-maf 0.001 --pbwt-modulo 0.025 --region $REG && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_shapeit5_$REG -y


done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt


