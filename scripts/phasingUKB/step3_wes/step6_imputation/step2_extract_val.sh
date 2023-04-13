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
	####################################
	#Setting 1
	ONAME="UKB_chr${CHR}.exome_array.WO_parents"
	#JOBID0=$(dx run app-swiss-army-knife -iin="/Phasing/PhasingWGS/step8_imputation/support/1k_samples_wbi.txt" -icmd="bcftools view /mnt/project/${PATHB}/chr${CHR}/${ONAME}.${METB2}.${EXTB} -S 1k_samples_wbi.txt -Ob -o wes_s1k_chr${CHR}.bcf --threads 4 && bcftools index -f wes_s1k_chr${CHR}.bcf --threads 4" --folder="/UKB_PHASING_EXOME_ARRAY/step6_imputation/validation" --instance-type mem1_ssd1_v2_x4 --priority low --name subs_${METB}_${CHR} -y | tail -n1 | cut -d" " -f3) 

	#JOBID0=$(dx run app-swiss-army-knife -icmd="bcftools view /mnt/project/UKB_PHASING_EXOME_ARRAY/step6_imputation/reference_panels/rp_s1k_${METS}_chr${CHR}.vcf.gz -G --threads 2 -Ob -o rp_s1k_${METS}_chr${CHR}_sites.bcf && bcftools index -f rp_s1k_${METS}_chr${CHR}_sites.bcf" --folder="/UKB_PHASING_EXOME_ARRAY/step6_imputation/reference_panels" --instance-type mem1_ssd1_v2_x2 --priority low --name sb_rp_${CHR} -y | tail -n1 | cut -d" " -f3) 

	JOBID0=$(dx run app-swiss-army-knife -icmd="bcftools view /mnt/project/UKB_PHASING_EXOME_ARRAY/step6_imputation/reference_panels/rp_s1k_${METS}_chr${CHR}.vcf.gz -S /mnt/project/data/ukb_wgs/ukb_cohort/wbic_samples_axiom_in_rp_wes.txt --threads 2 -G -Ob -o rp_wbic_${METS}_chr${CHR}_sites.bcf && bcftools index -f rp_wbic_${METS}_chr${CHR}_sites.bcf" --folder="/UKB_PHASING_EXOME_ARRAY/step6_imputation/reference_panels" --instance-type mem1_ssd1_v2_x2 --priority low --name sb_rp_af_${CHR} -y | tail -n1 | cut -d" " -f3) 

done
