#!/bin/bash

threads=16

for CHR in {1..22}; do
	ODIR=PhasingWGS_Official_release/step2_phase_rare/concat/chr${CHR}
	dx mkdir -p ${ODIR}	
	chunks=/mnt/project/PhasingWGS_Official_release/step2_phase_rare/chunks/UKB_chr${CHR}.chunk_*.shapeit5_rare.bcf
	OUT=UKB_chr${CHR}.full.shapeit5_rare.bcf
	
	dx run app-swiss-army-knife --folder="${ODIR}/" -icmd="ls -1v ${chunks} > concat_list_chr${CHR}.txt && bcftools concat -n -f concat_list_chr${CHR}.txt -o ${OUT} && bcftools index ${OUT} --threads ${threads} && rm concat_list_chr${CHR}.txt" --instance-type mem1_ssd1_v2_x16 --priority normal --name WGS_concat_chr${CHR} -y
	
done
