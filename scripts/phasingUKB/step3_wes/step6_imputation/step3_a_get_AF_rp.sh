#!/bin/bash

#Setting 1
PATHB="UKB_PHASING_EXOME_ARRAY/step5_phasing_beagle_common"
EXTB="common.vcf.gz.vcf.gz"
METB2="beagle"
METB="beagle5.4"

#Setting 2
PATHS="UKB_PHASING_EXOME_ARRAY/step3_phasing_rares"
EXTS="WO_parents.rares_0.001.bcf"
METS="shapeit5"

for CHR in {1..22}; do
	JOBID0=$(dx run app-swiss-army-knife -icmd="bcftools view /mnt/project/UKB_PHASING_EXOME_ARRAY/step6_imputation/reference_panels/rp_s1k_${METS}_chr${CHR}.vcf.gz -G -Ob -o rp_s1k_${METS}_chr${CHR}.bcf --threads 4 && bcftools index -f rp_s1k_${METS}_chr${CHR}.bcf --threads 4" --folder="/UKB_PHASING_EXOME_ARRAY/step6_imputation/reference_panels/sites" --instance-type mem1_ssd1_v2_x4 --priority low --name rp_pos_wes_${CHR} -y | tail -n1 | cut -d" " -f3) 

done
