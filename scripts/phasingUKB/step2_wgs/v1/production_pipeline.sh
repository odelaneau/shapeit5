#!/bin/bash

#INPUT FILES
BCF=/mnt/project/data/ukb_wgs/unphased/qc/ukb23352_c20_qc_v1.bcf
MAP=/mnt/project/data/shapeit_maps/chr20.b38.gmap.gz
PED=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/samples.families.caucasian.txt
CHK=/home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/v0/step2_splitchunks/chr20.size4Mb.txt

#BINARY
#DOCKER=shapeit5_$(git log -1 --format=%cd --date=short)\_$(git rev-parse --short HEAD)\.tar.gz
DOCKER=shapeit5.test.tar.gz

#STEP1: CONVERT TO SPARSE
COUT=ukb23352_c20_qc_v1.common.bcf
ROUT=ukb23352_c20_qc_v1.rare.bcf
dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step9_production/" -icmd="SHAPEIT5_scftools_static convert --input-plain $BCF --output-sparse $COUT $ROUT --region chr20 --thread 2 --maf 0.001" --tag bench1 --instance-type mem2_ssd1_v2_x2 --name bench1 --priority high -y

#STEP2: RUN COMMON VARIANT PHASING
while read LINE; do
	REG=$(echo $LINE | awk '{ print $3; }')
	BCF=ukb23352_c20_qc_v1.common.bcf
	OUT=ukb23352_c20_qc_v1.common.phased.$REG\.bcf
	LOG=ukb23352_c20_qc_v1.common.phased.$REG\.log
	#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step9_production/" -icmd="SHAPEIT5_phase_common_static --input $BCF --map $MAP --output $OUT --thread 32 --log $LOG --region $REG --pedigree $PED" --tag bench2 --instance-type mem3_ssd1_v2_x32 --priority high --name bench2_$REG -y
done < $CHK

#STEP3: LIGATE
rm ukb23352_c20_qc_v1.common.list
while read LINE; do
	REG=$(echo $LINE | awk '{ print $3; }')
	echo /mnt/project/Phasing/PhasingWGS/step9_production/ukb23352_c20_qc_v1.common.phased.$REG\.bcf >> ukb23352_c20_qc_v1.common.list
done < $CHK
dx cd /Phasing/PhasingWGS/step9_production/
dx rm ukb23352_c20_qc_v1.common.list
dx upload ukb23352_c20_qc_v1.common.list
TXT=ukb23352_c20_qc_v1.common.list
OUT=ukb23352_c20_qc_v1.common.phased.bcf
#dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step9_production/" -icmd="SHAPEIT5_ligate_static --input $TXT --output $OUT --thread 2 --index" --tag bench3 --instance-type mem2_ssd1_v2_x2 --priority high --name bench3 -y

#STEP4: RUN RARE VARIANT PHASING
while read LINE; do
	SREG=$(echo $LINE | awk '{ print $3; }')
	IREG=$(echo $LINE | awk '{ print $4; }')
	SPA=ukb23352_c20_qc_v1.rare.bcf
	SCA=ukb23352_c20_qc_v1.common.phased.bcf
	OUT=ukb23352_c20_qc_v1.rare.phased.$SREG\.bcf
	LOG=ukb23352_c20_qc_v1.rare.phased.$SREG\.log
	dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step9_production/" -icmd="SHAPEIT5_phase_rare_static --input-sparse $SPA --scaffold $SCA --map $MAP --pedigree $PED --output-sparse $OUT --log $LOG --scaffold-region $SREG --input-region $IREG --thread 32" --tag bench4 --instance-type mem3_ssd1_v2_x32 --priority normal --name bench4_$SREG -y
done < $CHK

#STEP5: MERGE SPARSE
rm ukb23352_c20_qc_v1.rare.list
while read LINE; do
	REG=$(echo $LINE | awk '{ print $3; }')
	echo /mnt/project/Phasing/PhasingWGS/step9_production/ukb23352_c20_qc_v1.rare.phased.$REG\.bcf >> ukb23352_c20_qc_v1.rare.list
done < $CHK
dx cd /Phasing/PhasingWGS/step9_production/
dx rm ukb23352_c20_qc_v1.rare.list
dx upload ukb23352_c20_qc_v1.rare.list
TXT=ukb23352_c20_qc_v1.rare.list
OUT=ukb23352_c20_qc_v1.rare.phased.bcf
dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step9_production/" -icmd="SHAPEIT5_scftools_static concat --input-sparse-list $TXT --output-sparse $OUT --thread 2" --tag bench5 --instance-type mem2_ssd1_v2_x2 --priority high --name bench5 -y

#STEP6: CONVERT TO BCF
CINP=ukb23352_c20_qc_v1.common.phased.bcf
RINP=ukb23352_c20_qc_v1.rare.phased.bcf
OUT=ukb23352_c20_qc_v1.phased.bcf
dx run app-swiss-army-knife -iimage_file="/docker/$DOCKER" --folder="/Phasing/PhasingWGS/step9_production/" -icmd="SHAPEIT5_scftools_static convert --input-sparse $CINP $RINP --output-plain $OUT --thread 2" --tag bench6 --instance-type mem2_ssd1_v2_x2 --priority high --name bench6 -y
