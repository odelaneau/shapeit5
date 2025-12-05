#!bin/bash


PREFIX=ukb200k
BIN=docker/shapeit5.ukb200k.tar.gz
threads=36
IDIR=Phasing/step2_ligate

ODIR=Phasing/step3_phase_rare/chunks
dx mkdir -p ${ODIR}

PED=/mnt/project/data/ukb_200k.ped


for CHR in X; do
	CHUNKS=../autosomes/chunks/phase_rare/chunks_chr${CHR}.txt
	BCF=/mnt/project/QC/chr${CHR}/ukb24304_chrX_without_PAR.Fmissing.bcf # QCeed data, unphased, full chromosome, full sample set.
	MAP=/mnt/project/data/shapeit_maps/chr${CHR}.b38.gmap.gz
	SCAF=/mnt/project/Phasing/step2_ligate/${PREFIX}_chr${CHR}.shapeit5_common_ligate.bcf

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


		OUT=${PREFIX}_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_rare.bcf
		LOG=${PREFIX}_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_rare.log
		TIM=${PREFIX}_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_rare.time


		echo "-------------------------------"
		echo "Chromosome: ${CHR}"
		echo "Chunk: ${CHUNK_NBR}"
		echo "Scaffold region: ${SCAFFOLD_REG}"
		echo "Input region: ${INPUT_REG}"
		echo "-------------------------------"


	dx run app-swiss-army-knife -iimage_file="/${BIN}" --folder="${ODIR}/" -icmd="/usr/bin/time -vo $TIM phase_rare_static --input $BCF --map $MAP --output $OUT --thread ${threads} --log $LOG --scaffold $SCAF --scaffold-region $SCAFFOLD_REG --input-region $INPUT_REG --pedigree ${PED} --haploids /mnt/project/data/chrX_males_to_keep_high_het_aneuploid_removed.txt && bcftools index -f $OUT --threads ${threads}" --instance-type mem1_ssd1_v2_x36 --priority normal --name WGS_shapeit5_rare_chr${CHR}_${CHUNK_NBR} -y
		
	done < ${CHUNKS}
done


