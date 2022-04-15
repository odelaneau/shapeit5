#!/bin/bash

#while read LINE; do
#	dx run app-swiss-army-knife --folder="/Phasing/PhasingWGS/step7_testPIR/" -icmd="cp $LINE ." --tag cp --instance-type mem2_ssd1_v2_x2 --priority low --name cp1 -y
#done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step7_testPIR/toyCRAMlist.txt

# WITH PIRS
while read LINE; do
	SREG=$(echo $LINE | awk '{ print $3; }')
	IREG=$(echo $LINE | awk '{ print $4; }')
	BCF=/mnt/project/Phasing/PhasingWGS/step2_splitchunks/benchmark_ukb23352_c20_qc_v1.$SREG\.bcf
	SCA=/mnt/project/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.default.depth8.bcf
	MAP=/mnt/project/data/shapeit_maps/chr20.b38.gmap.gz
	PIR=/mnt/project/data/ukb_wgs/support/cramfiles.paths.txt
	FAS=/mnt/project/data/reference_genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa
	
	OUT=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.phase2.pir.bcf
	LOG=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.phase2.pir.log
	TIM=benchmark_ukb23352_c20_qc_v1.$SREG\.shapeit5.phase2.pir.time
	
	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step6_runshapeit/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase2_static --input $BCF --scaffold $SCA --map $MAP --output $OUT --log $LOG --scaffold-region chr20:1000000-4000000 --input-region chr20:2000000-3000000 --bam-fasta $FAS --bam-list $PIR --thread 32 && bcftools index -f $OUT --threads 8" --tag benchWGS --tag shapeit5.2 --tag chr20 --instance-type mem3_ssd1_v2_x32 --priority low --name benchWGS_shapeit5.2_$SREG -y

done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.line1.txt



