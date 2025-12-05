#!/bin/bash

# Split the total number of phase_rare chunks to run several jobs in parrallel on the same machine.
n=8
mkdir -p chunks_rare/Splitted
rm chunks_rare/Splitted/*

cp chunks/phase_rare/chunks_chr*.txt chunks_rare/
cat chunks_rare/chunks_chr*.txt | cut -f -4 > chunks_rare/chunks.txt
split -d -a 3 -l ${n} chunks_rare/chunks.txt chunks_rare/Splitted/chunk_
N=$(($(ls chunks_rare/Splitted/chunk_* | wc -l)-1))


# phasing rare variants: 4 chunks_rare in parallel, 8 threads per job, 32 cpu total

PREFIX=ukb200k
BIN=docker/shapeit5.ukb200k.tar.gz
threads=8
IDIR=Phasing/step2_ligate
SDIR=Phasing/step3_phase_rare/support
ODIR=Phasing/step3_phase_rare/chunks
dx mkdir -p ${ODIR}
dx mkdir -p ${SDIR}
PED=ukb_200k.ped


SCR=script_phase_rare.sh
dx upload ${SCR} --path="${SDIR}/"



for SPL in $(seq 1 1 ${N}); do

	printf -v ID "%03d" $(echo $SPL)

	PARAM=chunks_rare/Splitted/chunk_${ID}
	dx rm "${SDIR}/$(basename ${PARAM})"
	dx upload ${PARAM} --path="${SDIR}/"
	PARAM=$(basename ${PARAM})
		
	dx run app-swiss-army-knife -iimage_file="/${BIN}" -iin="${SDIR}/${PARAM}" -iin="${SDIR}/${SCR}" -iin="data/${PED}" -icmd="cat ${PARAM} | xargs -P 4 -n4 bash ${SCR}" --instance-type mem3_ssd1_v2_x32 --folder="${ODIR}" --name shapeit5_rare_chunk_${SPL} --priority normal -y --brief --ignore-reuse


done
