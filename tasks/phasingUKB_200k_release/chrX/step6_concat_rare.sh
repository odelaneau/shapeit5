#!/bin/bash

PREFIX=ukb200k

threads=16

IDIR=Phasing/step3_phase_rare/chunks
ODIR=Phasing/step3_phase_rare

for CHR in X; do


	chunks=/mnt/project/${IDIR}/${PREFIX}_chr${CHR}.chunk_*.shapeit5_rare.bcf
	OUT=${PREFIX}_chr${CHR}.full.shapeit5_rare.bcf
	
	dx run app-swiss-army-knife --folder="${ODIR}/" -icmd="ls -1v ${chunks} > concat_list_chr${CHR}.txt && bcftools concat -n -f concat_list_chr${CHR}.txt -o ${OUT} && bcftools index ${OUT} --threads ${threads} && rm concat_list_chr${CHR}.txt" --instance-type mem1_ssd1_v2_x16 --priority normal --name WGS_concat_chr${CHR} -y
	
done
