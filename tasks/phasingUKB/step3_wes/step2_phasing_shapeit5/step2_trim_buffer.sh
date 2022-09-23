#!/bin/bash


for CHR in 14; do


	MAP=/mnt/project/data/shapeit_maps/chr${CHR}.b38.gmap.gz
#	CHUNKS=../step1_chunking/chunks.chr${CHR}.txt
	CHUNKS=../step1_chunking/robin_chunks_exome/chunks_chr${CHR}.txt
	CUT=0.001
	THREADS=16
	ODIR=UKB_PHASING_EXOME_ARRAY/step3_phasing_rares/chr${CHR}
	
	dx mkdir -p ${ODIR}


	while read LINE; do

			INPUT_REG=$(echo $LINE | awk '{ print $4; }')
			INPUT_REG_START=$(echo ${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f1)
			INPUT_REG_END=$(echo ${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f2)
			INPUT_REG_NAME=${CHR}_${INPUT_REG_START}_${INPUT_REG_END}
			INPUT_NBR=$(echo $LINE | awk '{ print $1; }')


			

			IN=/mnt/project/UKB_PHASING_EXOME_ARRAY/step3_phasing_rares/chr${CHR}/UKB_chr${CHR}.exome_array.${INPUT_REG_NAME}.shapeit5.WO_parents.rares_${CUT}.bcf
			
			OUT=UKB_chr${CHR}.exome_array.${INPUT_NBR}.shapeit5.WO_parents.rares_${CUT}.bcf
			echo $INPUT_NBR
			dx run app-swiss-army-knife --folder="./$ODIR" -icmd="bcftools view -r ${INPUT_REG} --threads ${THREADS} -Ob -o ${OUT} ${IN} && bcftools index ${OUT} --threads ${THREADS}" --tag shapeit5_rares --tag chr$CHR --tag $INPUT_REG_NAME --instance-type mem1_ssd1_v2_x16 --priority normal --name filter_$INPUT_REG_NAME -y
			
			
					
	done < $CHUNKS
done




