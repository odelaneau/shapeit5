#!/bin/bash

threads=16

for CHR in 14; do


	CUT=0.001
	ODIR=UKB_PHASING_EXOME_ARRAY/step3_phasing_rares

	dx mkdir -p ${ODIR}

	
	IN1=/mnt/project/${ODIR}/chr${CHR}/UKB_chr${CHR}.exome_array.1.shapeit5.WO_parents.rares_${CUT}.bcf
	IN2=/mnt/project/${ODIR}/chr${CHR}/UKB_chr${CHR}.exome_array.2.shapeit5.WO_parents.rares_${CUT}.bcf
	IN3=/mnt/project/${ODIR}/chr${CHR}/UKB_chr${CHR}.exome_array.3.shapeit5.WO_parents.rares_${CUT}.bcf
	
	
	OUT=UKB_chr${CHR}.exome_array.full.shapeit5.WO_parents.rares_${CUT}.bcf



	# CHR1-CHR4
	#dx run app-swiss-army-knife -icmd="bcftools concat --threads ${threads} --naive -Ob -o ${OUT} ${IN1} ${IN2} ${IN3} && bcftools index ${OUT} --threads ${threads}" --tag sites --instance-type mem1_ssd1_v2_x16 --folder="./${ODIR}/" --name merge_rares_chr${CHR} --priority normal -y

	# CHR5-CHR14
	dx run app-swiss-army-knife -icmd="bcftools concat --threads ${threads} --naive -Ob -o ${OUT} ${IN1} ${IN2} && bcftools index ${OUT} --threads ${threads}" --tag sites --instance-type mem1_ssd1_v2_x16 --folder="./${ODIR}/" --name merge_rares_chr${CHR} --priority normal -y
done




