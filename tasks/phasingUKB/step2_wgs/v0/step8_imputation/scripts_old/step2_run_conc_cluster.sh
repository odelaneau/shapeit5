#!/bin/bash

#FULL CHR
OFOLDER=/Phasing/PhasingWGS/step8_imputation/imputation/concordance_test

#for MET in beagle5.3 shapeit5; do
for MET in shapeit5; do
	#dx run app-swiss-army-knife -iimage_file="/docker/glimpse2_0.0.1.tar.gz" -icmd="echo \"chr20 /mnt/project/Phasing/PhasingWGS/step8_imputation/support/ukb_white_british_sites.bcf /mnt/project/Phasing/PhasingWGS/step8_imputation/validation/val_100_c20_qc_v1_all.bcf /mnt/project/Phasing/PhasingWGS/step8_imputation/imputation/imputed_100_rp_${MET}_chr20.vcf.gz\" > concordance.txt && GLIMPSE_concordance_v2.0.0 --gt-validation --minPROB 0 --minDP 0 --allele-counts --thread 2 --input concordance.txt --output imputed_100_rp_${MET} && rm -f concordance.txt" --tag concordance --instance-type mem1_ssd1_v2_x2 --priority normal --name ${MET}_conc --folder "${OFOLDER}" -y

	OFOLDER=/UKB_PHASING_WGS/step8_imputation/imputation/concordance_test
	dx run app-swiss-army-knife -iimage_file="/docker/glimpse2_0.0.1.tar.gz" -icmd="echo \"chr20 /mnt/project/Phasing/PhasingWGS/step8_imputation/support/ukb_white_british_sites.bcf /mnt/project/Phasing/PhasingWGS/step8_imputation/validation/val_100_c20_qc_v1_all.bcf /mnt/project/UKB_PHASING_WGS/step8_imputation/imputation/imputed_100_rp_${MET}_chr20.vcf.gz\" > concordance.txt && GLIMPSE_concordance_v2.0.0 --gt-validation --minPROB 0 --minDP 0 --allele-counts --thread 2 --input concordance.txt --output imputed_100_rp_${MET} && rm -f concordance.txt" --tag concordance --instance-type mem1_ssd1_v2_x2 --priority normal --name ${MET}_conc --folder "${OFOLDER}" -y
done
