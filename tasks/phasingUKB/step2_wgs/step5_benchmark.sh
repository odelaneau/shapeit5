#!/bin/bash

BCF=/mnt/project/Phasing/PhasingWGS/step1_preparedata/benchmark_ukb23352_c20_qc_v1.bcf
VAL=/mnt/project/Phasing/PhasingWGS/step1_preparedata/validation_ukb23352_c20_qc_v1.bcf
PED=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/samples.families.txt
FQA=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.SNPs.array.bcf
FQC=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.SNPs.common.bcf
FQR=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.SNPs.all.bcf
FQF=/mnt/project/Phasing/PhasingWGS/step1_preparedata/frequencies.ALL.all.bcf


#################################################################################################
#					GET FREQUENCY FILES					#
#################################################################################################

#FQCo=frequencies.SNPs.common.bcf
#FQRo=frequencies.SNPs.all.bcf
#FQFo=frequencies.ALL.all.bcf
#dx run app-swiss-army-knife --folder "/Phasing/PhasingWGS/step1_preparedata/" -icmd="bcftools view -G -Ob -o $FQCo -q 0.001:minor -v snps $BCF && bcftools index $FQCo" --tag getFRQ1 --tag benchWGS --instance-type mem2_ssd1_v2_x2 --name benchWGS_getFRQ1 --priority low -y
#dx run app-swiss-army-knife --folder "/Phasing/PhasingWGS/step1_preparedata/" -icmd="bcftools view -G -Ob -o $FQRo -v snps $BCF && bcftools index $FQRo" --tag getFRQ2 --tag benchWGS --instance-type mem2_ssd1_v2_x2 --name benchWGS_getFRQ2 --priority low -y
#dx run app-swiss-army-knife --folder "/Phasing/PhasingWGS/step1_preparedata/" -icmd="bcftools view -G -Ob -o $FQFo $BCF && bcftools index $FQFo" --tag getFRQ3 --tag benchWGS --instance-type mem2_ssd1_v2_x2 --name benchWGS_getFRQ3 --priority low -y

#################################################################################################
#				BENCHMARCH BEAGLE5						#
#################################################################################################

while read LINE; do
	iREG=$(echo $LINE | awk '{ print $3; }')
	oREG=$(echo $LINE | awk '{ print $4; }')
	#BEAGLE AT SNP ARRAY POSITIONS
	BGL=/mnt/project/Phasing/PhasingWGS/step3_runbeagle/benchmark_ukb23352_c20_qc_v1.$iREG\.beagle5.3.vcf.gz.vcf.gz
	OUT=benchmark_ukb23352_c20_qc_v1.$iREG\.beagle5.3.fqa
	LOG=benchmark_ukb23352_c20_qc_v1.$iREG\.beagle5.3.fqa.log
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step3_runbeagle/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BGL --frequency $FQA --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchBGL1_$iREG --tag chr20 --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchBGL1_$iREG -y
	
	#BEAGLE AT COMMON SNPs
	BGL=/mnt/project/Phasing/PhasingWGS/step3_runbeagle/benchmark_ukb23352_c20_qc_v1.$iREG\.beagle5.3.vcf.gz.vcf.gz
	OUT=benchmark_ukb23352_c20_qc_v1.$iREG\.beagle5.3.fqc
	LOG=benchmark_ukb23352_c20_qc_v1.$iREG\.beagle5.3.fqc.log
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step3_runbeagle/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BGL --frequency $FQC --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchBGL2_$iREG --tag chr20 --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchBGL2_$iREG -y
	
	#BEAGLE AT ALL SNPs
	BGL=/mnt/project/Phasing/PhasingWGS/step3_runbeagle/benchmark_ukb23352_c20_qc_v1.$iREG\.beagle5.3.vcf.gz.vcf.gz
	OUT=benchmark_ukb23352_c20_qc_v1.$iREG\.beagle5.3.fqr
	LOG=benchmark_ukb23352_c20_qc_v1.$iREG\.beagle5.3.fqr.log
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step3_runbeagle/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BGL --frequency $FQR --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchBGL3_$iREG --tag chr20 --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchBGL3_$iREG -y

	#BEAGLE AT ALL SITES
	BGL=/mnt/project/Phasing/PhasingWGS/step3_runbeagle/benchmark_ukb23352_c20_qc_v1.$iREG\.beagle5.3.vcf.gz.vcf.gz
	OUT=benchmark_ukb23352_c20_qc_v1.$iREG\.beagle5.3.fqf
	LOG=benchmark_ukb23352_c20_qc_v1.$iREG\.beagle5.3.fqf.log
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step3_runbeagle/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BGL --frequency $FQF --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchBGL4_$iREG --tag chr20 --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchBGL4_$iREG -y
done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt





#################################################################################################
#				BENCHMARCH SHAPEIT						#
#################################################################################################

while read LINE; do
	iREG=$(echo $LINE | awk '{ print $3; }')
	oREG=$(echo $LINE | awk '{ print $4; }')
	SHP=/mnt/project/Phasing/PhasingWGS/step6_runshapeit/benchmark_ukb23352_c20_qc_v1.$iREG\.shapeit5.phase2.bcf
	
	#SHAPEIT AT SNP ARRAY POSITIONS
	OUT=$(basename $SHP)\.fqa
	LOG=$(basename $SHP)\.fqa.log
	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step6_runshapeit/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $SHP --frequency $FQA --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP1_$iREG --tag chr20 --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchSHP1_$iREG -y
	
	#SHAPEIT AT COMMON SNPs
	OUT=$(basename $SHP)\.fqc
	LOG=$(basename $SHP)\.fqc.log
	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step6_runshapeit/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $SHP --frequency $FQC --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP2_$iREG --tag chr20 --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchSHP2_$iREG -y
	
	#SHAPEIT AT ALL SNPs
	OUT=$(basename $SHP)\.fqr
	LOG=$(basename $SHP)\.fqr.log
	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step6_runshapeit/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $SHP --frequency $FQR --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP3_$iREG --tag chr20 --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchSHP3_$iREG -y

	#SHAPEIT AT ALL SITES
	OUT=$(basename $SHP)\.fqf
	LOG=$(basename $SHP)\.fqf.log
	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step6_runshapeit/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $SHP --frequency $FQF --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP4_$iREG --tag chr20 --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchSHP4_$iREG -y
done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt







#################################################################################################
#				BENCHMARCH SHAPEIT5 SCAFFOLDS					#
#################################################################################################


while read LINE; do
	iREG=$(echo $LINE | awk '{ print $3; }')
	oREG=$(echo $LINE | awk '{ print $4; }')

	########################################################## DEPTH 4 MODULO 0.1cM ############################################################################

	#SHAPEIT AT SNP ARRAY POSITIONS
	SHP=/mnt/project/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.$iREG\.shapeit5.default.bcf
	OUT=$(basename $SHP .bcf)\.fqa
	LOG=$(basename $SHP .bcf)\.fqa.log
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $SHP --frequency $FQA --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP1_$iREG --tag chr20 --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchSHP1_$iREG -y
	
	#SHAPEIT AT COMMON SNPs
	SHP=/mnt/project/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.$iREG\.shapeit5.default.bcf
	OUT=$(basename $SHP .bcf)\.fqc
	LOG=$(basename $SHP .bcf)\.fqc.log
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $SHP --frequency $FQC --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP2_$iREG --tag chr20 --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchSHP2_$iREG -y
	
	
	#SHAPEIT AT SNP ARRAY POSITIONS
	SHP=/mnt/project/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.$iREG\.shapeit5.scaffold.bcf
	OUT=$(basename $SHP .bcf)\.fqa
	LOG=$(basename $SHP .bcf)\.fqa.log
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $SHP --frequency $FQA --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP3_$iREG --tag chr20 --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchSHP3_$iREG -y
	
	#SHAPEIT AT COMMON SNPs
	SHP=/mnt/project/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.$iREG\.shapeit5.scaffold.bcf
	OUT=$(basename $SHP .bcf)\.fqc
	LOG=$(basename $SHP .bcf)\.fqc.log
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $SHP --frequency $FQC --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP4_$iREG --tag chr20 --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchSHP4_$iREG -y
	
	########################################################## DEPTH 8 MODULO 0.1cM ############################################################################
	
	#SHAPEIT AT SNP ARRAY POSITIONS
	SHP=/mnt/project/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.$iREG\.shapeit5.default.depth8.bcf
	OUT=$(basename $SHP .bcf)\.fqa
	LOG=$(basename $SHP .bcf)\.fqa.log
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $SHP --frequency $FQA --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP1_$iREG --tag chr20 --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchSHP1_$iREG -y
	
	#SHAPEIT AT COMMON SNPs
	SHP=/mnt/project/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.$iREG\.shapeit5.default.depth8.bcf
	OUT=$(basename $SHP .bcf)\.fqc
	LOG=$(basename $SHP .bcf)\.fqc.log
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $SHP --frequency $FQC --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP2_$iREG --tag chr20 --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchSHP2_$iREG -y
	
	
	#SHAPEIT AT SNP ARRAY POSITIONS
	SHP=/mnt/project/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.$iREG\.shapeit5.scaffold.depth8.bcf
	OUT=$(basename $SHP .bcf)\.fqa
	LOG=$(basename $SHP .bcf)\.fqa.log
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $SHP --frequency $FQA --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP3_$iREG --tag chr20 --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchSHP3_$iREG -y
	
	#SHAPEIT AT COMMON SNPs
	SHP=/mnt/project/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.$iREG\.shapeit5.scaffold.depth8.bcf
	OUT=$(basename $SHP .bcf)\.fqc
	LOG=$(basename $SHP .bcf)\.fqc.log
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $SHP --frequency $FQC --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP4_$iREG --tag chr20 --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchSHP4_$iREG -y

	########################################################## DEPTH 4 MODULO 0.05cM ############################################################################
	
	#SHAPEIT AT SNP ARRAY POSITIONS
	SHP=/mnt/project/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.$iREG\.shapeit5.default.modulo5.bcf
	OUT=$(basename $SHP .bcf)\.fqa
	LOG=$(basename $SHP .bcf)\.fqa.log
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $SHP --frequency $FQA --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP1_$iREG --tag chr20 --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchSHP1_$iREG -y
	
	#SHAPEIT AT COMMON SNPs
	SHP=/mnt/project/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.$iREG\.shapeit5.default.modulo5.bcf
	OUT=$(basename $SHP .bcf)\.fqc
	LOG=$(basename $SHP .bcf)\.fqc.log
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $SHP --frequency $FQC --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP2_$iREG --tag chr20 --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchSHP2_$iREG -y
	
	
	#SHAPEIT AT SNP ARRAY POSITIONS
	SHP=/mnt/project/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.$iREG\.shapeit5.scaffold.modulo5.bcf
	OUT=$(basename $SHP .bcf)\.fqa
	LOG=$(basename $SHP .bcf)\.fqa.log
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $SHP --frequency $FQA --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP3_$iREG --tag chr20 --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchSHP3_$iREG -y
	
	#SHAPEIT AT COMMON SNPs
	SHP=/mnt/project/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.$iREG\.shapeit5.scaffold.modulo5.bcf
	OUT=$(basename $SHP .bcf)\.fqc
	LOG=$(basename $SHP .bcf)\.fqc.log
	#dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_0.0.1.tar.gz" --folder="/Phasing/PhasingWGS/step4_runshapeit/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $SHP --frequency $FQC --pedigree $PED --region $oREG --output $OUT --log $LOG --thread 2" --tag benchWGS --tag switchSHP4_$iREG --tag chr20 --instance-type mem3_ssd1_v2_x2 --priority low --name benchWGS_switchSHP4_$iREG -y

	
done < /home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt



