#!/bin/bash

# step0. Create a docker folder of the RAP and upload the provided shapeit5 docker image

dx mkdir -p docker/
dx upload ../../docker/shapeit5_v1.0.0.tar.gz --path="docker/"
BIN=docker/shapeit5_v1.0.0.tar.gz

# step0. Create output directory.
ODIR=PhasingWGS_Official_release/step1_phase_common/chunks
dx mkdir -p ${ODIR}

# step1. Upload map files
cp ../../maps/genetic_maps.b38.tar.gz ./
tar xvzf genetic_maps.b38.tar.gz
dx mkdir -p data/shapeit_maps/
dx upload *.b38.gmap.gz --path="data/shapeit_maps/"


# step2. Get phasing chunks
cp ../../chunks/chunks.b38.tar.gz ./
tar -xvzf chunks.b38.tar.gz


# step3. phasing common variants (MAF>0.001)
for CHR in {1..22}; do

	CHUNKS=chunks/phase_common/chunks_chr${CHR}.txt
	BCF=/mnt/project/PhasingWGS_Official_release/step0_qc/ukb23352_c${CHR}_qc_v1.bcf # QCeed data from step0, unphased, full chromosome, full sample set.
	MAP=/mnt/project/data/shapeit_maps/chr${CHR}.b38.gmap.gz

	while read LINE; do
		REG=$(echo $LINE | awk '{ print $3; }')
		CHUNK_NBR=$(echo $LINE | awk '{ print $1; }')
		OUT=UKB_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_common.bcf
		LOG=UKB_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_common.log
		TIM=UKB_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_common.time
		
		dx run app-swiss-army-knife -iimage_file="/${BIN}" --folder="${ODIR}/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase_common_static_v1.0.0 --input $BCF --map $MAP --output $OUT --thread 72 --log $LOG --filter-maf 0.001 --region $REG && bcftools index -f $OUT --threads 72" --instance-type mem1_ssd1_v2_x72 --priority normal --name WGS_shapeit5_common_chr${CHR}_${CHUNK_NBR} -y
		
	done < ${CHUNKS}
done



