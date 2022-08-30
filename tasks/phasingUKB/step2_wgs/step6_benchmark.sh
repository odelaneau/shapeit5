#!/bin/bash

BCF=/mnt/project/Phasing/PhasingWGS/step1_preparedata/benchmark_ukb23352_c20_qc_v1.bcf
VAL=/mnt/project/Phasing/PhasingWGS/step1_preparedata/validation_ukb23352_c20_qc_v1.bcf
PED=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/samples.families.txt

DOCKER=shapeit5_$(git log -1 --format=%cd --date=short)\_$(git rev-parse --short HEAD)\.tar.gz

#####################################################################################################
#				BENCHMARCH BEAGLE5							#
#####################################################################################################
for N in 2000 5000 10000 20000 50000 100000 147754; do
	while read LINE; do
		iREG=$(echo $LINE | awk '{ print $3; }')
		oREG=$(echo $LINE | awk '{ print $4; }')
		VCF=/mnt/project/Phasing/PhasingWGS/step3_runbeagle/N$N/benchmark_ukb23352_c20_qc_v1.subset.N$N\.$iREG\.vcf.gz
		
		#SNP ARRAY POSITIONS
		OUT=$(basename $VCF .vcf.gz)\.fqa
		LOG=$(basename $VCF .vcf.gz)\.fqa.log
		FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.SNPs.array.bcf
		#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step3_runbeagle/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $VCF --frequency $FRQ --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchBGL1 --tag $iREG --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchBGL1_$iREG -y
		
		#COMMON SNPS
		OUT=$(basename $VCF .vcf.gz)\.fqc
		LOG=$(basename $VCF .vcf.gz)\.fqc.log
		FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.SNPs.common.bcf
		#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step3_runbeagle/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $VCF --frequency $FRQ --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchBGL2 --tag $iREG --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchBGL2_$iREG -y
		
		#RARE SNPS
		OUT=$(basename $VCF .vcf.gz)\.fqr
		LOG=$(basename $VCF .vcf.gz)\.fqr.log
		FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.SNPs.all.bcf
		#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step3_runbeagle/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $VCF --frequency $FRQ --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchBGL3 --tag $iREG --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchBGL3_$iREG -y
		
		#ALL VARIANTS
		OUT=$(basename $VCF .vcf.gz)\.fqf
		LOG=$(basename $VCF .vcf.gz)\.fqf.log
		FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.ALL.all.bcf
		#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step3_runbeagle/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $VCF --frequency $FRQ --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchBGL4 --tag $iREG --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchBGL4_$iREG -y
	
	done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt
done


#################################################################################################
#				BENCHMARCH SHAPEIT5						#
#################################################################################################

#for N in 2000 5000 10000 20000 50000 100000 147754; do
for N in 147754; do
	while read LINE; do
		#for T in default scaffold; do
		for T in default; do
			iREG=$(echo $LINE | awk '{ print $3; }')
			oREG=$(echo $LINE | awk '{ print $4; }')
			BCF=/mnt/project/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/benchmark_ukb23352_c20_qc_v1.subset.N$N\.$iREG\.shapeit5.$T\.bcf
			
			#SNP ARRAY POSITIONS
			OUT=$(basename $BCF .bcf)\.fqa
			LOG=$(basename $BCF .bcf)\.fqa.log
			FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.SNPs.array.bcf
			#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BCF --frequency $FRQ --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP1 --tag $iREG --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchSHP1_$iREG -y
			
			#COMMON SNPS
			OUT=$(basename $BCF .bcf)\.fqc
			LOG=$(basename $BCF .bcf)\.fqc.log
			FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.SNPs.common.bcf
			#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BCF --frequency $FRQ --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP2 --tag $iREG --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchSHP2_$iREG -y
			
			#RARE SNPS
			OUT=$(basename $BCF .bcf)\.fqr
			LOG=$(basename $BCF .bcf)\.fqr.log
			FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.SNPs.all.bcf
			#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BCF --frequency $FRQ --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP3 --tag $iREG --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchSHP3_$iREG -y
			
			#ALL VARIANTS
			OUT=$(basename $BCF .bcf)\.fqf
			LOG=$(basename $BCF .bcf)\.fqf.log
			FRQ=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.subset.N$N\.ALL.all.bcf
			dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step5_runshapeit_phase2/N$N/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BCF --frequency $FRQ --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP4 --tag $iREG --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchSHP4_$iREG -y
		
		done
	done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt
done

