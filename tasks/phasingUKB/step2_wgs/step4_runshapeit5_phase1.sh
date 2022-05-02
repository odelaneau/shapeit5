#!/bin/bash

MAP=/mnt/project/data/shapeit_maps/chr20.b38.gmap.gz
SCA=/mnt/project/Phasing/PhasingSNParray/step5_benchmark/benchmark_c20_b0_v2.b38.sorted.N480853.phased.bcf

for N in 2000 5000 10000 20000 50000 100000 147754; do

	while read LINE; do
		
		REG=$(echo $LINE | awk '{ print $3; }')
		BCF=/mnt/project/Phasing/PhasingWGS/step2_splitchunks/benchmark_ukb23352_c20_qc_v1.subset.N$N\.$REG\.bcf
		
		OUT=benchmark_ukb23352_c20_qc_v1.subset.N$N\.$REG\.shapeit5.default.bcf
		LOG=benchmark_ukb23352_c20_qc_v1.subset.N$N\.$REG\.shapeit5.default.log
		TIM=benchmark_ukb23352_c20_qc_v1.subset.N$N\.$REG\.shapeit5.default.time
		
		dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit_phase1/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase1_static --input $BCF --map $MAP --output $OUT --thread 32 --log $LOG --filter-maf 0.001 --region $REG && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_shapeit5_$REG -y
		
		OUT=benchmark_ukb23352_c20_qc_v1.subset.N$N\.$REG\.shapeit5.scaffold.bcf
		LOG=benchmark_ukb23352_c20_qc_v1.subset.N$N\.$REG\.shapeit5.scaffold.log
		TIM=benchmark_ukb23352_c20_qc_v1.subset.N$N\.$REG\.shapeit5.scaffold.time
		
		dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit_phase1/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase1_static --input $BCF --scaffold $SCA --map $MAP --output $OUT --thread 32 --log $LOG --filter-maf 0.001 --region $REG && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_shapeit5_$REG -y
	
	done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt
done

