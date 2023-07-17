#!bin/bash

PREFIX=ukb200k
PED=ukb_200k.ped
threads=8

CHR=$(echo $2 | sed 's/chr//')

BCF=/mnt/project/QC/chr${CHR}/ukb24304_chr${CHR}.qceed.bcf # QCeed data, unphased, full chromosome, full sample set.
MAP=/mnt/project/data/shapeit_maps/chr${CHR}.b38.gmap.gz
SCAF=/mnt/project/Phasing/step2_ligate/${PREFIX}_chr${CHR}.shapeit5_common_ligate.bcf

CHUNK_NBR=$1
SCAFFOLD_REG=$3
SCAFFOLD_REG_START=$(echo ${SCAFFOLD_REG} | cut -d":" -f 2 | cut -d"-" -f1)
SCAFFOLD_REG_END=$(echo ${SCAFFOLD_REG} | cut -d":" -f 2 | cut -d"-" -f2)
SCAFFOLD_REG_NAME=${CHR}_${SCAFFOLD_REG_START}_${SCAFFOLD_REG_END}
			
INPUT_REG=$4
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


/usr/bin/time -vo $TIM phase_rare_static --input $BCF --map $MAP --output $OUT --thread ${threads} --log $LOG --scaffold $SCAF --scaffold-region $SCAFFOLD_REG --input-region $INPUT_REG --pedigree ${PED} && bcftools index -f $OUT --threads ${threads}



