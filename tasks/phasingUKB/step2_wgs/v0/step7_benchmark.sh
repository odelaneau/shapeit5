#!/bin/bash

BCF=/mnt/project/Phasing/PhasingWGS/step1_preparedata/benchmark_ukb23352_c20_qc_v1.bcf
VAL=/mnt/project/Phasing/PhasingWGS/step1_preparedata/validation_ukb23352_c20_qc_v1.bcf
PEDALL=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/samples.families.caucasian.txt
PEDDUO=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/samples.duos.caucasian.txt
PEDTRI=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/samples.trios.caucasian.txt

DOCKER=shapeit5.test.tar.gz

#####################################################################################################
#				BENCHMARCH BEAGLE5		CHUNKS 														#
#####################################################################################################

for N in 2000 5000 10000 20000 50000 100000 147754; do
	VCF=/mnt/project/Phasing/PhasingWGS/step3_runbeagle/N$N/benchmark_ukb23352_c20_qc_v1.subset.N$N\.fullchr.vcf.gz
	
	while read LINE; do
		REG=$(echo $LINE | awk '{ print $4; }')
		
		#SNP ARRAY POSITIONS
		OUT=$(basename $VCF .vcf.gz)\.$REG\.fqa
		LOG=$(basename $VCF .vcf.gz)\.$REG\.fqa.log
		FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.SNPs.array.bcf
		#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step3_runbeagle/N$N/chunks/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $VCF --frequency $FRQ --pedigree $PEDALL --region $REG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchBGL1 --instance-type mem2_ssd1_v2_x2 --priority normal --name benchWGS_switchBGL1 -y
		
		#COMMON SNPS
		OUT=$(basename $VCF .vcf.gz)\.$REG\.fqc
		LOG=$(basename $VCF .vcf.gz)\.$REG\.fqc.log
		FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.SNPs.common.bcf
		#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step3_runbeagle/N$N/chunks/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $VCF --frequency $FRQ --pedigree $PEDALL --region $REG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchBGL2 --instance-type mem3_ssd1_v2_x2 --priority normal --name benchWGS_switchBGL2 -y
			
		#ALL VARIANTS
		OUT=$(basename $VCF .vcf.gz)\.$REG\.fqf
		LOG=$(basename $VCF .vcf.gz)\.$REG\.fqf.log
		FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.ALL.all.bcf
		#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step3_runbeagle/N$N/chunks/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $VCF --frequency $FRQ --pedigree $PEDALL --region $REG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchBGL3 --instance-type mem3_ssd1_v2_x2 --priority normal --name benchWGS_switchBGL3 -y
			
		#ALL VARIANTS - DUOS
		OUT=$(basename $VCF .vcf.gz)\.$REG\.fqd
		LOG=$(basename $VCF .vcf.gz)\.$REG\.fqd.log
		FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.ALL.all.bcf
		#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step3_runbeagle/N$N/chunks/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $VCF --frequency $FRQ --pedigree $PEDDUO --region $REG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchBGL4 --instance-type mem3_ssd1_v2_x2 --priority normal --name benchWGS_switchBGL4 -y
			
		#ALL VARIANTS - TRIOS
		OUT=$(basename $VCF .vcf.gz)\.$REG\.fqt
		LOG=$(basename $VCF .vcf.gz)\.$REG\.fqt.log
		FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.ALL.all.bcf
		#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step3_runbeagle/N$N/chunks/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $VCF --frequency $FRQ --pedigree $PEDTRI --region $REG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchBGL5 --instance-type mem3_ssd1_v2_x2 --priority normal --name benchWGS_switchBGL5 -y
	done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/v0/step2_splitchunks/chr20.size4Mb.txt
done

#################################################################################################
#				BENCHMARCH SHAPEIT5	 CHUNKS													#
#################################################################################################

for N in 2000 5000 10000 20000 50000 100000 147754; do
	BCF=/mnt/project/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/benchmark_ukb23352_c20_qc_v1.subset.N$N\.fullchr.shapeit5.ligated.bcf
	
	while read LINE; do
		REG=$(echo $LINE | awk '{ print $4; }')	
		
		#SNP ARRAY POSITIONS
		OUT=$(basename $BCF .bcf)\.$REG\.fqa
		LOG=$(basename $BCF .bcf)\.$REG\.fqa.log
		FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.SNPs.array.bcf
		#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/chunks/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BCF --frequency $FRQ --pedigree $PEDALL --region $REG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP1 --instance-type mem2_ssd1_v2_x2 --priority normal --name benchWGS_switchSHP1 -y
				
		#COMMON SNPS
		OUT=$(basename $BCF .bcf)\.$REG\.fqc
		LOG=$(basename $BCF .bcf)\.$REG\.fqc.log
		FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.SNPs.common.bcf
		#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/chunks/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BCF --frequency $FRQ --pedigree $PEDALL --region $REG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP2 --instance-type mem3_ssd1_v2_x2 --priority normal --name benchWGS_switchSHP2 -y
		
		#ALL VARIANTS
		OUT=$(basename $BCF .bcf)\.$REG\.fqf
		LOG=$(basename $BCF .bcf)\.$REG\.fqf.log
		FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.ALL.all.bcf
		#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/chunks/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BCF --frequency $FRQ --pedigree $PEDALL --region $REG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP3 --instance-type mem3_ssd1_v2_x2 --priority normal --name benchWGS_switchSHP3 -y
		
		#ALL VARIANTS - DUOS
		OUT=$(basename $BCF .bcf)\.$REG\.fqd
		LOG=$(basename $BCF .bcf)\.$REG\.fqd.log
		FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.ALL.all.bcf
		#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/chunks/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BCF --frequency $FRQ --pedigree $PEDDUO --region $REG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP4 --instance-type mem3_ssd1_v2_x2 --priority normal --name benchWGS_switchSHP4 -y
		
		#ALL VARIANTS - TRIOS
		OUT=$(basename $BCF .bcf)\.$REG\.fqt
		LOG=$(basename $BCF .bcf)\.$REG\.fqt.log
		FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.ALL.all.bcf
		#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/chunks/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BCF --frequency $FRQ --pedigree $PEDTRI --region $REG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP5 --instance-type mem3_ssd1_v2_x2 --priority normal --name benchWGS_switchSHP5 -y
	done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/v0/step2_splitchunks/chr20.size4Mb.txt		

done




for N in 147754; do
	
	while read LINE; do
		SREG=$(echo $LINE | awk '{ print $3; }')
		IREG=$(echo $LINE | awk '{ print $4; }')
		BCF=/mnt/project/Phasing/PhasingWGS/step5_runshapeit_phase2/N147754.test/benchmark_ukb23352_c20_qc_v1.subset.N147754.$SREG\.shapeit5.default.bcf
		#BCF=/mnt/project/Phasing/PhasingWGS/step9_production/ukb23352_c20_qc_v1.common.phased.bcf
		#BCF=/mnt/project/phasing_rare/rephased/N147754/benchmark_ukb23352_c20_qc_v1.subset.N147754.fullchr.shapeit5.ligated.bcf_rephased.bcf 
				         
		#ALL VARIANTS
		OUT=$(basename $BCF .bcf)\.$IREG\.default.fqf
		LOG=$(basename $BCF .bcf)\.$IREG\.default.fqf.log
		FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.ALL.all.bcf
		dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N147754.test/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BCF --frequency $FRQ --pedigree $PEDALL --region $IREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP3 --instance-type mem3_ssd1_v2_x2 --priority normal --name benchWGS_switchSHP3 -y
	done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/v0/step2_splitchunks/chr20.size4Mb.txt		

done
