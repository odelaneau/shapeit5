#!/bin/bash


# To phase the chromosome X, we designed a specific --haploid option that takes as input the list of haploid individuals (i.e males). To determine a confident set of haploid males, we used several QC steps that are described in the phasing report (https://docs.google.com/document/d/1EJmh-JcR8HBvu3rjBtREw_50kDIW-sOW2EcC6zweuTc/edit?usp=sharing). The haploid individuals are stored in the file /mnt/project/data/chrX_males_to_keep_high_het_aneuploid_removed.txt.


PREFIX=ukb200k

BIN=docker/shapeit5.ukb200k.tar.gz


# step0. Create output directory.
ODIR=Phasing/step1_phase_common/chunks
dx mkdir -p ${ODIR}

# step1. Upload map files
#cp ../Git_repository/shapeit5/maps/genetic_maps.b38.tar.gz ./
#tar xvzf genetic_maps.b38.tar.gz
#dx mkdir -p data/shapeit_maps/
#dx upload *.b38.gmap.gz --path="data/shapeit_maps/"


PED=/mnt/project/data/ukb_200k.ped


# step2. phasing common variants (MAF>0.001)
for CHR in X; do


	CHUNKS=../autosomes/chunks/phase_common/chunks_chr${CHR}.txt
	BCF=/mnt/project/QC/chr${CHR}/ukb24304_chrX_without_PAR.bcf # QCeed data, unphased, full chromosome, full sample set.
	MAP=/mnt/project/data/shapeit_maps/chr${CHR}.b38.gmap.gz

	while read LINE; do
		REG=$(echo $LINE | awk '{ print $3; }')
		CHUNK_NBR=$(echo $LINE | awk '{ print $1; }')
		OUT=${PREFIX}_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_common.bcf
		LOG=${PREFIX}_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_common.log
		TIM=${PREFIX}_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_common.time
		
		dx run app-swiss-army-knife -iimage_file="/${BIN}" --folder="${ODIR}/" -icmd="/usr/bin/time -vo $TIM phase_common_static --input $BCF --map $MAP --output $OUT --thread 72 --log $LOG --filter-maf 0.001 --region $REG --pedigree ${PED} --haploid /mnt/project/data/chrX_males_to_keep_high_het_aneuploid_removed.txt && bcftools index -f $OUT --threads 72" --instance-type mem1_ssd1_v2_x72 --priority normal --name ${PREFIX}_shapeit5_common_chr${CHR}_${CHUNK_NBR} -y --brief --ignore-reuse
		
	done < ${CHUNKS}
done






