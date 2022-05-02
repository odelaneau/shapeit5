#!/bin/bash


MAP=/mnt/project/data/shapeit_maps/chr20.b38.gmap.gz
SCA=/mnt/project/Phasing/PhasingSNParray/step5_benchmark/benchmark_c20_b0_v2.b38.sorted.N480853.phased.bcf

for N in 2000 5000 10000 20000 50000 100000 147754; do

	while read LINE; do
		
		SREG=$(echo $LINE | awk '{ print $3; }')
		IREG=$(echo $LINE | awk '{ print $4; }')
		BCF=/mnt/project/Phasing/PhasingWGS/step2_splitchunks/benchmark_ukb23352_c20_qc_v1.subset.N$N\.$REG\.bcf

		#NOT SCAFFOLDED		
		SCA=/mnt/project/Phasing/PhasingWGS/step4_runshapeit_phase1/benchmark_ukb23352_c20_qc_v1.subset.N$N\.$REG\.shapeit5.default.bcf
		OUT=benchmark_ukb23352_c20_qc_v1.subset.N$N\.$REG\.shapeit5.default
		
		dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/" -icmd="/usr/bin/time -vo $OUT\.time SHAPEIT5_phase2_static --input $BCF --scaffold $SCA --map $MAP --output $OUT\.bcf --log $OUT\.log --scaffold-region $SREG --input-region $IREG --thread 32 && bcftools index -f $OUT\.bcf --threads 8" --tag benchWGS --tag shapeit5 --tag phase2 --tag $SREG --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_shapeit5_phase2a_$SREG -y
		
		#SCAFFOLDED
		SCA=/mnt/project/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.subset.N$N\.$REG\.shapeit5.scaffold.bcf
		OUT=benchmark_ukb23352_c20_qc_v1.subset.N$N\.$REG\.shapeit5.scaffold
		dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/" -icmd="/usr/bin/time -vo $OUT\.time SHAPEIT5_phase2_static --input $BCF --scaffold $SCA --map $MAP --output $OUT\.bcf --log $OUT\.log --scaffold-region $SREG --input-region $IREG --thread 32 && bcftools index -f $OUT\.bcf --threads 8" --tag benchWGS --tag shapeit5 --tag phase2 --tag $SREG --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_shapeit5_phase2b_$SREG -y
		
	done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt

done
