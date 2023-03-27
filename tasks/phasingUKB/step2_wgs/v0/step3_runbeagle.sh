#!/bin/bash


BGL=/mnt/project/docker/beagle.19Apr22.7c0.jar
MAP=/mnt/project/data/plink_maps/plink.prefix.chr20.GRCh38.map

#BY CHUNKS
for N in 2000 5000 10000 20000 50000 100000 147754; do

	while read LINE; do
		
		REG=$(echo $LINE | awk '{ print $3; }')
		VCF=/mnt/project/Phasing/PhasingWGS/step2_splitchunks/benchmark_ukb23352_c20_qc_v1.subset.N$N\.$REG\.vcf.gz
		OUT=$(basename $VCF .vcf.gz)
		TIM=$(basename $VCF .vcf.gz)\.time
		
		#dx run app-swiss-army-knife --folder="/Phasing/PhasingWGS/step3_runbeagle/N$N/" -icmd="/usr/bin/time -vo $TIM java -Xmx256G -jar $BGL gt=$VCF map=$MAP out=$OUT nthreads=32 chrom=chr20 && bcftools index -f $OUT\.vcf.gz --threads 8" --tag benchWGS --tag beagle5.4 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority normal --name benchWGS_beagle5.4_$N\_$REG -y
	
	done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt
done

#WHOLE CHROMOSOME
for N in 2000 5000 10000 20000 50000 100000; do
	VCF=/mnt/project/Phasing/PhasingWGS/step2_splitchunks/benchmark_ukb23352_c20_qc_v1.subset.N$N\.vcf.gz
	OUT=$(basename $VCF .vcf.gz)\.fullchr
	TIM=$(basename $VCF .vcf.gz)\.fullchr.time
	
	dx run app-swiss-army-knife --folder="/Phasing/PhasingWGS/step3_runbeagle/N$N/" -icmd="/usr/bin/time -vo $TIM java -Xmx460G -jar $BGL gt=$VCF map=$MAP out=$OUT nthreads=64 window=20 chrom=chr20 && bcftools index -f $OUT\.vcf.gz --threads 8" --tag benchWGS --tag beagle5.4 --tag chr20 --instance-type mem3_ssd1_v2_x64 --priority high --name benchWGS_beagle5.4 -y
done
