#!/bin/bash


BIN=docker/shapeit5_v1.0.0.tar.gz
threads=36

# step0. Create output directory.
ODIR=PhasingWGS_Official_release/step2_phase_rare/chunks
dx mkdir -p ${ODIR}

# step1. phasing rare variants (MAF<0.001)
for CHR in {1..22}; do
	CHUNKS=chunks/phase_rare/chunks_chr${CHR}.txt
	BCF=/mnt/project/PhasingWGS_Official_release/step0_qc/ukb23352_c${CHR}_qc_v1.bcf # QCeed data from step0, unphased, full chromosome, full sample set.
	MAP=/mnt/project/data/shapeit_maps/chr${CHR}.b38.gmap.gz
	SCAF=/mnt/project/PhasingWGS_Official_release/step1_phase_common/UKB_chr${CHR}.shapeit5_common_ligate.bcf	
	
	while read LINE; do

		CHUNK_NBR=$(echo $LINE | awk '{ print $1; }')
		
		SCAFFOLD_REG=$(echo $LINE | awk '{ print $3; }')
		SCAFFOLD_REG_START=$(echo ${SCAFFOLD_REG} | cut -d":" -f 2 | cut -d"-" -f1)
		SCAFFOLD_REG_END=$(echo ${SCAFFOLD_REG} | cut -d":" -f 2 | cut -d"-" -f2)
		SCAFFOLD_REG_NAME=${CHR}_${SCAFFOLD_REG_START}_${SCAFFOLD_REG_END}
			
		INPUT_REG=$(echo $LINE | awk '{ print $4; }')
		INPUT_REG_START=$(echo ${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f1)
		INPUT_REG_END=$(echo ${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f2)
		INPUT_REG_NAME=${CHR}_${INPUT_REG_START}_${INPUT_REG_END}

		OUT=UKB_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_rare.bcf
		LOG=UKB_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_rare.log
		TIM=UKB_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_rare.time
		
		dx run app-swiss-army-knife -iimage_file="/${BIN}" --folder="${ODIR}/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase_rare_static_v1.0.0 --input-plain $BCF --map $MAP --output $OUT --thread ${threads} --log $LOG --scaffold $SCAF --scaffold-region $SCAFFOLD_REG --input-region $INPUT_REG && bcftools index -f $OUT --threads ${threads}" --instance-type mem3_ssd1_v2_x48 --priority normal --name WGS_shapeit5_rare_chr${CHR}_${CHUNK_NBR} -y
		
	done < ${CHUNKS}
done



