#!/bin/bash

for CHR in {1..22}; do
	VCF=/mnt/project/data/SNParray_hg19_phased_shapeit4/UKB.chr${CHR}.snps_QC_filter.phased_shapeit4.bcf
	OUT=UKB.chr${CHR}.snps_QC_filter.phased_shapeit4.indices.repbwt
	TIM=UKB.chr${CHR}.snps_QC_filter.phased_shapeit4.indices.repbwt.time
	dx run app-swiss-army-knife --folder "/data/SNParray_hg19_phased_shapeit4/repbwt_indices/" -iimage_file="/docker/repbwt_abe107e_20230119.tar.gz" -icmd="/usr/bin/time -vo $TIM repbwt -i ${VCF} -s ${OUT}" --instance-type mem3_ssd1_v2_x8 --priority normal --name repbwt_chr${CHR} -y
done

