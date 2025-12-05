#!/bin/bash

CHR=20
VCF=/mnt/project/Phasing/PhasingSNParray/step4_liftover/full_c20_b0_v2.b38.sorted.vcf.gz

#GET VALIDATION
SUB=/mnt/project/Phasing/PhasingSNParray/step5_benchmark/validation.samples.txt
OUT=validation_c20_b0_v2.b38.sorted.vcf.gz
#dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step5_benchmark/" -icmd="bcftools view -Oz -o $OUT -S $SUB --force-samples $VCF && bcftools index $OUT" --tag subset --tag benchmark --instance-type mem2_ssd1_v2_x2 --name benchmark_subset --priority normal -y

#GET DOWNSAMPLES
for N in 5000 10000 20000 50000 100000 200000 300000 400000 480853; do
	OUT=benchmark_c20_b0_v2.b38.sorted.N$N\.vcf.gz
	SUB=/mnt/project/Phasing/PhasingSNParray/step5_benchmark/N$N\.txt
	#dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step5_benchmark/" -icmd="bcftools view -Oz -o $OUT -S $SUB $VCF && bcftools index $OUT" --tag subset --tag benchmark --instance-type mem2_ssd1_v2_x2 --name benchmark_subset --priority normal -y
done

#PHASING SHAPEIT
MAP=/mnt/project/data/shapeit_maps/chr20.b38.gmap.gz
for N in 5000 10000 20000 50000 100000 200000 300000 400000 480853; do
	INP=/mnt/project/Phasing/PhasingSNParray/step5_benchmark/benchmark_c20_b0_v2.b38.sorted.N$N\.vcf.gz
	LOG=benchmark_c20_b0_v2.b38.sorted.N$N\.phased.log
	OUT=benchmark_c20_b0_v2.b38.sorted.N$N\.phased.bcf
	TIM=benchmark_c20_b0_v2.b38.sorted.N$N\.phased.time
	#dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step5_benchmark/" -iimage_file="/docker/shapeit5_0.0.1.tar.gz" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase1_static --input $INP --map $MAP --output $OUT --region chr20 --log $LOG --thread 16 && bcftools index $OUT --threads 16" --tag phasing --tag benchmark --instance-type mem2_ssd1_v2_x16 --priority normal --name benchmark_phasing -y
done

#PHASING BEAGLE
BGL=/mnt/project/docker/beagle.19Apr22.7c0.jar
MAP=/mnt/project/data/plink_maps/plink.prefix.chr20.GRCh38.map
for N in 5000 10000 20000 50000 100000 200000 300000 400000 480853; do
	INP=/mnt/project/Phasing/PhasingSNParray/step5_benchmark/benchmark_c20_b0_v2.b38.sorted.N$N\.vcf.gz
	OUT=benchmark_c20_b0_v2.b38.sorted.N$N\.beagle
	#dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step5_benchmark/" -icmd="/usr/bin/time -vo $OUT\.time java -Xmx64G -jar $BGL gt=$INP map=$MAP out=$OUT nthreads=16 chrom=chr20 && bcftools index -f $OUT\.vcf.gz --threads 16" --tag phasing --tag benchmark --tag beagle5.4 --instance-type mem2_ssd1_v2_x16 --priority normal --name benchmark_phasing_beagle5.4 -y
done

#BENCHMARK SHAPEIT
VAL=/mnt/project/Phasing/PhasingSNParray/step5_benchmark/validation_c20_b0_v2.b38.sorted.vcf.gz
PED=/mnt/project/Phasing/PhasingSNParray/step5_benchmark/families.fam
for N in 5000 10000 20000 50000 100000 200000 300000 400000 480853; do
	INP=/mnt/project/Phasing/PhasingSNParray/step5_benchmark/benchmark_c20_b0_v2.b38.sorted.N$N\.phased.bcf
	FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.SNPs.array.bcf
	OUT=benchmark_c20_b0_v2.b38.sorted.N$N\.phased.fqa
	LOG=benchmark_c20_b0_v2.b38.sorted.N$N\.phased.fqa.log
	
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingSNParray/step5_benchmark/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $INP --frequency $FRQ --pedigree $PED --region chr20 --output $OUT --log $LOG --thread 2" --tag benchSNP --tag switchSHP --tag $N --instance-type mem3_ssd1_v2_x2 --priority low --name benchSNP_switchSHP_$N -y
done

#BENCHMARK SHAPEIT USING SEQUENCING DATA
VAL=/mnt/project/Phasing/PhasingWGS/step1_preparedata/validation_ukb23352_c20_qc_v1.bcf
PED=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/samples.families.txt
for N in 5000 10000 20000 50000 100000 200000 300000 400000 480853; do
	INP=/mnt/project/Phasing/PhasingSNParray/step5_benchmark/benchmark_c20_b0_v2.b38.sorted.N$N\.phased.bcf
	FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.SNPs.array.bcf
	OUT=benchmark_c20_b0_v2.b38.sorted.N$N\.phased.wgs.fqa
	LOG=benchmark_c20_b0_v2.b38.sorted.N$N\.phased.wgs.fqa.log
	
	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingSNParray/step5_benchmark/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $INP --frequency $FRQ --pedigree $PED --region chr20 --output $OUT --log $LOG --thread 2 --dupid" --tag benchSNP --tag switchSHP --tag $N --instance-type mem3_ssd1_v2_x2 --priority low --name benchSNP_switchSHP_$N -y
done

#BENCHMARK BEAGLE
VAL=/mnt/project/Phasing/PhasingSNParray/step5_benchmark/validation_c20_b0_v2.b38.sorted.vcf.gz
PED=/mnt/project/Phasing/PhasingSNParray/step5_benchmark/families.fam
for N in 5000 10000 20000 50000 100000 200000 300000 400000 480853; do
	INP=/mnt/project/Phasing/PhasingSNParray/step5_benchmark/benchmark_c20_b0_v2.b38.sorted.N$N\.beagle.vcf.gz
	FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.SNPs.array.bcf
	OUT=benchmark_c20_b0_v2.b38.sorted.N$N\.beagle.fqa
	LOG=benchmark_c20_b0_v2.b38.sorted.N$N\.beagle.fqa.log
	
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingSNParray/step5_benchmark/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $INP --frequency $FRQ --pedigree $PED --region chr20 --output $OUT --log $LOG --thread 2" --tag benchSNP --tag switchBGL --tag $N --instance-type mem3_ssd1_v2_x2 --priority low --name benchSNP_switchBGL_$N -y
done

#RETRIEVING SHAPEIT SERs
#dx cd PhasingUKB:/Phasing/PhasingSNParray/step5_benchmark
#for N in 5000 10000 20000 50000 100000 200000 300000 400000 480853; do 
#	dx cat benchmark_c20_b0_v2.b38.sorted.N$N\.phased.fqa.sample.mendel.txt.gz | zcat | grep -v "\-1" | awk '{ print $1; }' > samples.txt;
#	vT=$(dx cat benchmark_c20_b0_v2.b38.sorted.N$N\.phased.fqa.sample.switch.txt.gz | zcat | grep -wf samples.txt  | awk 'BEGIN {e=0;t=0; } { t+=$3; e+=$2; } END { print e*100/t; }');
#	vD=$(dx cat benchmark_c20_b0_v2.b38.sorted.N$N\.phased.fqa.sample.switch.txt.gz | zcat | awk 'BEGIN {e=0;t=0; } { t+=$3; e+=$2; } END { print e*100/t; }');
#	rm samples.txt;
#	echo $N $vT $vD;
#done

