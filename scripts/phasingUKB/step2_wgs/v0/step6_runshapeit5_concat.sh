#!/bin/bash

for N in 2000 5000 10000 20000 50000 100000 147754; do
	
	#Produce list of files
	rm benchmark_ukb23352_c20_qc_v1.subset.N$N\.shapeit5.default.list
	while read LINE; do
		REG=$(echo $LINE | awk '{ print $3; }')
		echo /mnt/project/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/benchmark_ukb23352_c20_qc_v1.subset.N$N\.$REG\.shapeit5.default.bcf >> benchmark_ukb23352_c20_qc_v1.subset.N$N\.shapeit5.default.list
	done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt
	dx cd /Phasing/PhasingWGS/step5_runshapeit_phase2/N$N
	dx rm benchmark_ukb23352_c20_qc_v1.subset.N$N\.shapeit5.default.list
	dx upload benchmark_ukb23352_c20_qc_v1.subset.N$N\.shapeit5.default.list
	rm benchmark_ukb23352_c20_qc_v1.subset.N$N\.shapeit5.default.list
	
	#Run job
	TXT=/mnt/project/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/benchmark_ukb23352_c20_qc_v1.subset.N$N\.shapeit5.default.list
	OUT=benchmark_ukb23352_c20_qc_v1.subset.N$N\.fullchr.shapeit5.default.bcf
	
	dx run app-swiss-army-knife --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/" -icmd="bcftools concat -n -f $TXT -Ob -o $OUT --threads 2 && bcftools index -f $OUT --threads 2" --tag benchWGS --tag shapeit5 --tag phase2 --instance-type mem2_ssd1_v2_x2 --priority normal --name benchWGS_concat -y

done

