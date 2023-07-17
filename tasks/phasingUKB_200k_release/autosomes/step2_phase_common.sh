#!/bin/bash

CHR_START=$1
CHR_END=$2


PREFIX=ukb200k

BIN=docker/shapeit5.ukb200k.tar.gz # if not done yet, upload the SHAPEIT5 docker image onto the RAP.


# step0. Create output directory.
ODIR=Phasing/step1_phase_common/chunks
dx mkdir -p ${ODIR}

# step1. Upload map files
S5_DIR=../Git_repository/shapeit5 ## specify your own path for the SHAPEIT5 github
cp ${S5_DIR}/maps/genetic_maps.b38.tar.gz ./
tar xvzf genetic_maps.b38.tar.gz
dx mkdir -p data/shapeit_maps/
dx upload *.b38.gmap.gz --path="data/shapeit_maps/"


PED=/mnt/project/data/ukb_200k.ped # the pedifree file contains three columns: col1=offspring; col2=father; col3=mother. NA can be specified for unknown parent. Parent-offspring trios and duos are determined from the UK Biobank relatedness file (provided as part of the official released - computed using the KING software). We used 0.1717<Kinship<0.3535 & IBS0<0.0012.


# step2. phasing common variants (MAF>0.001)
for CHR in $(seq $CHR_START 1 $CHR_END); do


	CHUNKS=chunks/phase_common/chunks_chr${CHR}.txt
	BCF=/mnt/project/QC/chr${CHR}/ukb24304_chr${CHR}.qceed.bcf # QCeed data, unphased, full chromosome, full sample set.
	MAP=/mnt/project/data/shapeit_maps/chr${CHR}.b38.gmap.gz

	while read LINE; do
		REG=$(echo $LINE | awk '{ print $3; }')
		CHUNK_NBR=$(echo $LINE | awk '{ print $1; }')
		OUT=${PREFIX}_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_common.bcf
		LOG=${PREFIX}_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_common.log
		TIM=${PREFIX}_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_common.time
		
		dx run app-swiss-army-knife -iimage_file="/${BIN}" --folder="${ODIR}/" -icmd="/usr/bin/time -vo $TIM phase_common_static --input $BCF --map $MAP --output $OUT --thread 72 --log $LOG --filter-maf 0.001 --region $REG --pedigree ${PED} && bcftools index -f $OUT --threads 72" --instance-type mem1_ssd1_v2_x72 --priority normal --name ${PREFIX}_shapeit5_common_chr${CHR}_${CHUNK_NBR} -y --brief --ignore-reuse
		
	done < ${CHUNKS}
done



