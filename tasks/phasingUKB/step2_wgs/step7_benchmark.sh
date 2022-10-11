#!/bin/bash

BCF=/mnt/project/Phasing/PhasingWGS/step1_preparedata/benchmark_ukb23352_c20_qc_v1.bcf
VAL=/mnt/project/Phasing/PhasingWGS/step1_preparedata/validation_ukb23352_c20_qc_v1.bcf
PEDALL=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/samples.families.caucasian.txt
PEDDUO=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/samples.duos.txt
PEDTRI=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/samples.trios.txt

DOCKER=shapeit5_$(git log -1 --format=%cd --date=short)\_$(git rev-parse --short HEAD)\.tar.gz

#####################################################################################################
#				BENCHMARCH BEAGLE5		CHUNKS 														#
#####################################################################################################

#for N in 2000 5000 10000 20000 50000 100000 147754; do
for N in 147754; do
	VCF=/mnt/project/Phasing/PhasingWGS/step3_runbeagle/N$N/benchmark_ukb23352_c20_qc_v1.subset.N$N\.fullchr.vcf.gz
	
	while read LINE; do
		REG=$(echo $LINE | awk '{ print $4; }')
		
		#SNP ARRAY POSITIONS
		OUT=$(basename $VCF .vcf.gz)\.$REG\.fqa
		LOG=$(basename $VCF .vcf.gz)\.$REG\.fqa.log
		FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.SNPs.array.bcf
		dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step3_runbeagle/$N/chunks/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $VCF --frequency $FRQ --pedigree $PEDALL --region $REG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchBGL1 --instance-type mem3_ssd1_v2_x2 --priority normal --name benchWGS_switchBGL1 -y
		
		#COMMON SNPS
		OUT=$(basename $VCF .vcf.gz)\.$REG\.fqc
		LOG=$(basename $VCF .vcf.gz)\.$REG\.fqc.log
		FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.SNPs.common.bcf
		dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step3_runbeagle/$N/chunks/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $VCF --frequency $FRQ --pedigree $PEDALL --region $REG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchBGL2 --instance-type mem3_ssd1_v2_x2 --priority normal --name benchWGS_switchBGL2 -y
			
		#ALL VARIANTS
		OUT=$(basename $VCF .vcf.gz)\.$REG\.fqf
		LOG=$(basename $VCF .vcf.gz)\.$REG\.fqf.log
		FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.ALL.all.bcf
		dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step3_runbeagle/$N/chunks/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $VCF --frequency $FRQ --pedigree $PEDALL --region $REG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchBGL3 --instance-type mem3_ssd1_v2_x2 --priority normal --name benchWGS_switchBGL3 -y
	
	done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt
done

#####################################################################################################
#				BENCHMARCH BEAGLE5		FULLCHR														#
#####################################################################################################

#for N in 2000 5000 10000 20000 50000 100000 147754; do
for N in 147754; do
	VCF=/mnt/project/Phasing/PhasingWGS/step3_runbeagle/N$N/benchmark_ukb23352_c20_qc_v1.subset.N$N\.fullchr.vcf.gz
	
	#SNP ARRAY POSITIONS
	OUT=$(basename $VCF .vcf.gz)\.fqa
	LOG=$(basename $VCF .vcf.gz)\.fqa.log
	FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.SNPs.array.bcf
	#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step3_runbeagle/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $VCF --frequency $FRQ --pedigree $PEDALL --region chr20 --output $OUT --log $LOG --thread 8" --tag benchWGS --tag switchBGL1 --instance-type mem3_ssd1_v2_x8 --priority normal --name benchWGS_switchBGL1 -y
		
	#COMMON SNPS
	OUT=$(basename $VCF .vcf.gz)\.fqc
	LOG=$(basename $VCF .vcf.gz)\.fqc.log
	FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.SNPs.common.bcf
	#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step3_runbeagle/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $VCF --frequency $FRQ --pedigree $PEDALL --region chr20 --output $OUT --log $LOG --thread 8" --tag benchWGS --tag switchBGL2 --instance-type mem3_ssd1_v2_x8 --priority normal --name benchWGS_switchBGL2 -y
		
	#ALL VARIANTS
	OUT=$(basename $VCF .vcf.gz)\.fqf
	LOG=$(basename $VCF .vcf.gz)\.fqf.log
	FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.ALL.all.bcf
	#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step3_runbeagle/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $VCF --frequency $FRQ --pedigree $PEDALL --region chr20 --output $OUT --log $LOG --thread 8" --tag benchWGS --tag switchBGL3 --instance-type mem3_ssd1_v2_x8 --priority normal --name benchWGS_switchBGL3 -y

done


#################################################################################################
#				BENCHMARCH SHAPEIT5	 FULL CHR													#
#################################################################################################

for N in 2000 5000 10000 20000 50000 100000 147754; do
	BCF=/mnt/project/Phasing/PhasingWGS/step5_runshapeit_phase2/benchmark_ukb23352_c20_qc_v1.subset.N$N\.fullchr.shapeit5.default.bcf	
	
	#SNP ARRAY POSITIONS
	OUT=$(basename $BCF .bcf)\.fqa
	LOG=$(basename $BCF .bcf)\.fqa.log
	FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.SNPs.array.bcf
	#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BCF --frequency $FRQ --pedigree $PEDALL --region chr20 --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP1 --instance-type mem3_ssd1_v2_x2 --priority normal --name benchWGS_switchSHP1 -y
			
	#COMMON SNPS
	OUT=$(basename $BCF .bcf)\.fqc
	LOG=$(basename $BCF .bcf)\.fqc.log
	FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.SNPs.common.bcf
	#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BCF --frequency $FRQ --pedigree $PEDALL --region chr20 --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP2 --instance-type mem3_ssd1_v2_x2 --priority normal --name benchWGS_switchSHP2 -y
		
	#ALL VARIANTS
	OUT=$(basename $BCF .bcf)\.fqf
	LOG=$(basename $BCF .bcf)\.fqf.log
	FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.ALL.all.bcf
	#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BCF --frequency $FRQ --pedigree $PEDALL --region chr20 --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP3 --instance-type mem3_ssd1_v2_x2 --priority normal --name benchWGS_switchSHP3 -y
		
	#SINGLETONS : TRIOS ONLY
	OUT=$(basename $BCF .bcf)\.fqt
	LOG=$(basename $BCF .bcf)\.fqt.log
	FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.ALL.all.bcf
	#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BCF --frequency $FRQ --pedigree $PEDTRI --region chr20 --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP4 --instance-type mem3_ssd1_v2_x2 --priority normal --name benchWGS_switchSHP4 -y
		
	#SINGLETONS : DUOS ONLY + MENDEL
	OUT=$(basename $BCF .bcf)\.fqf.duos0
	LOG=$(basename $BCF .bcf)\.fqf.duos0.log
	FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.ALL.all.bcf
	#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BCF --frequency $FRQ --pedigree $PEDDUO --region chr20 --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP5 --instance-type mem3_ssd1_v2_x2 --priority normal --name benchWGS_switchSHP5 -y
		
	#SINGLETONS : DUOS ONLY + --singleton
	OUT=$(basename $BCF .bcf)\.fqf.duos1
	LOG=$(basename $BCF .bcf)\.fqf.duos1.log
	FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.ALL.all.bcf
	#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BCF --frequency $FRQ --pedigree $PEDDUO --region chr20 --output $OUT --singleton --log $LOG --thread 2" --tag benchWGS --tag switchSHP6 --instance-type mem3_ssd1_v2_x2 --priority normal --name benchWGS_switchSHP6 -y

done

#################################################################################################
#				BENCHMARCH SHAPEIT5	 CHUNKS													#
#################################################################################################

for N in 2000 5000 10000 20000 50000 100000 147754; do
	BCF=/mnt/project/Phasing/PhasingWGS/step5_runshapeit_phase2/benchmark_ukb23352_c20_qc_v1.subset.N$N\.fullchr.shapeit5.default.bcf
	
	while read LINE; do
		REG=$(echo $LINE | awk '{ print $4; }')	
	
		#SNP ARRAY POSITIONS
		OUT=$(basename $BCF .bcf)\.$REG\.fqa
		LOG=$(basename $BCF .bcf)\.$REG\.fqa.log
		FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.SNPs.array.bcf
		#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BCF --frequency $FRQ --pedigree $PEDALL --region chr20 --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP1 --instance-type mem3_ssd1_v2_x2 --priority normal --name benchWGS_switchSHP1 -y
				
		#COMMON SNPS
		OUT=$(basename $BCF .bcf)\.$REG\.fqc
		LOG=$(basename $BCF .bcf)\.$REG\.fqc.log
		FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.SNPs.common.bcf
		#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BCF --frequency $FRQ --pedigree $PEDALL --region chr20 --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP2 --instance-type mem3_ssd1_v2_x2 --priority normal --name benchWGS_switchSHP2 -y
		
		#ALL VARIANTS
		OUT=$(basename $BCF .bcf)\.$REG\.fqf
		LOG=$(basename $BCF .bcf)\.$REG\.fqf.log
		FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.ALL.all.bcf
		#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BCF --frequency $FRQ --pedigree $PEDALL --region chr20 --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP3 --instance-type mem3_ssd1_v2_x2 --priority normal --name benchWGS_switchSHP3 -y
	done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt		

done

