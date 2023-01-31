#!/bin/bash

CHR=20
for CNK in {0..12}; do
	VCF=/mnt/project/PhasingWGS_Official_release/step2_phase_rare/chunks/UKB_chr${CHR}.chunk_${CNK}.shapeit5_rare.bcf
	OUT=UKB_chr${CHR}.chunk_${CNK}.shapeit5_rare.indices.repbwt
	TIM=UKB_chr${CHR}.chunk_${CNK}.shapeit5_rare.indices.repbwt.time
	dx run app-swiss-army-knife --folder "/data/SNParray_hg19_phased_shapeit4/repbwt_indices_WGS_chr20/" -iimage_file="/docker/repbwt_abe107e_20230119.tar.gz" -icmd="/usr/bin/time -vo $TIM repbwt -i ${VCF} -s ${OUT}" --instance-type mem3_ssd1_v2_x8 --priority normal --name repbwt_chr${CHR}_wgs_${CNK} -y
done

