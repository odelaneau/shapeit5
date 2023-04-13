#!/bin/bash


for CHR in {1..22}; do


	MAP=/mnt/project/data/shapeit_maps/chr${CHR}.b38.gmap.gz
	CHUNKS=../step1_chunking/chunk_wes/chunks_chr${CHR}.txt
	CUT=0.001
	THREADS=32
	ODIR=UKB_PHASING_EXOME_ARRAY/step3_phasing_rares/chr${CHR}
	
	dx mkdir -p ${ODIR}


	while read LINE; do

			SCAFFOLD_REG=$(echo $LINE | awk '{ print $3; }')
			SCAFFOLD_REG_START=$(echo ${SCAFFOLD_REG} | cut -d":" -f 2 | cut -d"-" -f1)
			SCAFFOLD_REG_END=$(echo ${SCAFFOLD_REG} | cut -d":" -f 2 | cut -d"-" -f2)
			SCAFFOLD_REG_NAME=${CHR}_${SCAFFOLD_REG_START}_${SCAFFOLD_REG_END}

			INPUT_REG=$(echo $LINE | awk '{ print $4; }')
			INPUT_REG_START=$(echo ${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f1)
			INPUT_REG_END=$(echo ${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f2)
			INPUT_REG_NAME=${CHR}_${INPUT_REG_START}_${INPUT_REG_END}

			BCF=/mnt/project/UKB_PHASING_EXOME_ARRAY/step0_merge/chr${CHR}/UKB.chr${CHR}.exome_array.WO_parents.sorted.bcf

			echo ${SCAFFOLD_REG}; echo ${INPUT_REG}
			
			SCAFFOLD=/mnt/project/UKB_PHASING_EXOME_ARRAY/step2_phasing_common/chr${CHR}/UKB_chr${CHR}.exome_array.WO_parents.shapeit5.common_${CUT}.bcf
			OUT=UKB_chr${CHR}.exome_array.${INPUT_REG_NAME}.shapeit5.WO_parents.rares_${CUT}.bcf
			LOG=UKB_chr${CHR}.exome_array.${INPUT_REG_NAME}.shapeit5.WO_parents.rares_${CUT}.log
			TIM=UKB_chr${CHR}.exome_array.${INPUT_REG_NAME}.shapeit5.WO_parents.rares_${CUT}.time

			dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_2022-08-30_12bfbd7.tar.gz" --folder="./$ODIR" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase_rare_static --input $BCF --scaffold $SCAFFOLD --map $MAP --output $OUT --log $LOG --scaffold-region $SCAFFOLD_REG --input-region $INPUT_REG --thread $THREADS && bcftools index -f $OUT --threads $THREADS" --tag shapeit5_rares --tag chr$CHR --tag $INPUT_REG_NAME --instance-type mem2_ssd1_v2_x32 --priority normal --name shapeit5_rares_$INPUT_REG_NAME -y
			
			
					
	done < $CHUNKS
done




