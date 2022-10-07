#!/bin/bash

OFOLDER=/Phasing/PhasingWGS/step8_imputation/target_snp_array

for CHR in {1..22}; do
	dx run app-swiss-army-knife -icmd="bcftools view -S /mnt/project/Phasing/PhasingWGS/step8_imputation/support/1k_samples_imp_wbi.txt /mnt/project/Phasing/PhasingSNParray/step4_liftover/full_c${CHR}_b0_v2.b38.sorted.vcf.gz --threads 2 -Oz -o t_axiom_1k_c${CHR}_b0_v2.b38.vcf.gz && bcftools reheader -s /mnt/project/Phasing/PhasingWGS/step8_imputation/support/rename_1k_samples_wbi.txt t_axiom_1k_c${CHR}_b0_v2.b38.vcf.gz -o axiom_1k_c${CHR}_b0_v2.b38.vcf.gz && bcftools index -f axiom_1k_c${CHR}_b0_v2.b38.vcf.gz --threads 2 && rm -f t_axiom_1k_c${CHR}_b0_v2.b38.vcf.gz" --instance-type mem1_ssd1_v2_x2 --priority low --name s8_s0_subs_array_chr${CHR} --folder "${OFOLDER}" -y
done

