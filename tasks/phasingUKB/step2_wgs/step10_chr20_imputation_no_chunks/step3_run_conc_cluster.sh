#!/bin/bash

#FULL CHR
OFOLDER=/Phasing/PhasingWGS/step8_imputation/concordance_test/binning2

METB="beagle5.4"
METS="shapeit5"

########### OLD
	#dx run app-swiss-army-knife -iimage_file="/docker/glimpse2_20220825_502047a.tar.gz" -icmd="echo \"chr20 /mnt/project/data/ukb_wgs/unphased/step1_wbic_noval/ukb23352_c20_qc_v1_wbic_no10kgreml_sites.bcf /mnt/project/data/ukb_wgs/phased_shapeit/panels/validation_greml/wgs_ukb23352_c20_qc_v1_all.bcf /mnt/project/Phasing/PhasingWGS/step8_imputation/ligated/imputed_1k_rp_${METB}_chr${CHR}.bcf\" > concordance.txt && GLIMPSE_concordance_v2.0.0 --gt-val --allele-counts --thread 2 --input concordance.txt --output imputed_1k_rp_${METB}_chr${CHR} && rm -f concordance.txt" --instance-type mem1_ssd1_v2_x2 --priority low --name ${METB}_conc_wgs --folder "${OFOLDER}" -y

	#dx run app-swiss-army-knife -iimage_file="/docker/glimpse2_20220825_502047a.tar.gz" -icmd="echo \"chr20 /mnt/project/data/ukb_wgs/unphased/step1_wbic_noval/ukb23352_c20_qc_v1_wbic_no10kgreml_sites.bcf /mnt/project/data/ukb_wgs/phased_shapeit/panels/validation_greml/wgs_ukb23352_c20_qc_v1_all.bcf /mnt/project/Phasing/PhasingWGS/step8_imputation/ligated/imputed_1k_rp_${METS}_chr${CHR}.bcf\" > concordance.txt && GLIMPSE_concordance_v2.0.0 --gt-val --allele-counts --thread 2 --input concordance.txt --output imputed_1k_rp_${METS}_chr${CHR} && rm -f concordance.txt" --instance-type mem1_ssd1_v2_x2 --priority low --name ${METS}_conc_wgs --folder "${OFOLDER}" -y

	#dx run app-swiss-army-knife -iimage_file="/docker/glimpse2_20221002_73decec-dirty.tar.gz" -icmd="echo \"chr20 /mnt/project/Phasing/PhasingWGS/step8_imputation/reference_panels/sites/rp_s1k_shapeit5_chr20.bcf /mnt/project/data/ukb_wgs/phased_shapeit/panels/validation_greml/wgs_ukb23352_c20_qc_v1_all.bcf /mnt/project/Phasing/PhasingWGS/step8_imputation/ligated/imputed_1k_rp_${METB}_chr${CHR}.bcf\" > concordance.txt && GLIMPSE_concordance_v2.0.0 --gt-val --allele-counts --thread 2 --input concordance.txt --output imputed_1k_rp_${METB}_chr${CHR}_rp_af && rm -f concordance.txt" --instance-type mem1_ssd1_v2_x2 --priority low --name ${METB}_conc_wgs --folder "${OFOLDER}" -y

	#dx run app-swiss-army-knife -iimage_file="/docker/glimpse2_20221002_73decec-dirty.tar.gz" -icmd="echo \"chr20 /mnt/project/Phasing/PhasingWGS/step8_imputation/reference_panels/sites/rp_s1k_shapeit5_chr20.bcf /mnt/project/data/ukb_wgs/phased_shapeit/panels/validation_greml/wgs_ukb23352_c20_qc_v1_all.bcf /mnt/project/Phasing/PhasingWGS/step8_imputation/ligated/imputed_1k_rp_${METS}_chr${CHR}.bcf\" > concordance.txt && GLIMPSE_concordance_v2.0.0 --gt-val --allele-counts --thread 2 --input concordance.txt --output imputed_1k_rp_${METS}_chr${CHR}_rp_af && rm -f concordance.txt" --instance-type mem1_ssd1_v2_x2 --priority low --name ${METS}_conc_wgs --folder "/Phasing/PhasingWGS/step8_imputation/concordance_test/" -y

	#dx run app-swiss-army-knife -iimage_file="/docker/glimpse2_20221005_73decec-dirty.tar.gz" -icmd="echo \"chr20 /mnt/project/data/ukb_wgs/unphased/step1_wbic_noval/ukb23352_c20_qc_v1_wbic_no10kgreml_sites.bcf /mnt/project/data/ukb_wgs/phased_shapeit/panels/validation_greml/wgs_ukb23352_c20_qc_v1_all.bcf /mnt/project/Phasing/PhasingWGS/step8_imputation/ligated/imputed_1k_rp_${METB}_chr${CHR}.bcf\" > concordance.txt && GLIMPSE_concordance_v2.0.0 --gt-val --ac-bins 1 5 10 20 50 100 200 500 1000  2000 5000 10000  20000 50000 100000 146754 --input concordance.txt --output imputed_1k_rp_${METB}_chr${CHR}_wbi_af_binning2 && rm -f concordance.txt" --instance-type mem1_ssd1_v2_x2 --priority low --name ${METB}_conc_wgs --folder "${OFOLDER}" -y

	#dx run app-swiss-army-knife -iimage_file="/docker/glimpse2_20221005_73decec-dirty.tar.gz" -icmd="echo \"chr20 /mnt/project/data/ukb_wgs/unphased/step1_wbic_noval/ukb23352_c20_qc_v1_wbic_no10kgreml_sites.bcf /mnt/project/data/ukb_wgs/phased_shapeit/panels/validation_greml/wgs_ukb23352_c20_qc_v1_all.bcf /mnt/project/Phasing/PhasingWGS/step8_imputation/ligated/imputed_1k_rp_${METS}_chr${CHR}.bcf\" > concordance.txt && GLIMPSE_concordance_v2.0.0 --gt-val --ac-bins 1 5 10 20 50 100 200 500 1000  2000 5000 10000  20000 50000 100000 146754 --thread 2 --input concordance.txt --output imputed_1k_rp_${METS}_chr${CHR}_wbi_af_binning2 && rm -f concordance.txt" --instance-type mem1_ssd1_v2_x2 --priority low --name ${METS}_conc_wgs --folder "${OFOLDER}" -y


##########BINNING 2 RP - NEW, for the paper
	#dx run app-swiss-army-knife -iimage_file="/docker/glimpse2_20221005_73decec-dirty.tar.gz" -icmd="echo \"chr20 /mnt/project/Phasing/PhasingWGS/step8_imputation/reference_panels/sites/rp_s1k_shapeit5_chr20.bcf /mnt/project/data/ukb_wgs/phased_shapeit/panels/validation_greml/wgs_ukb23352_c20_qc_v1_all.bcf /mnt/project/Phasing/PhasingWGS/step8_imputation/ligated/imputed_1k_rp_${METB}_chr${CHR}.bcf\" > concordance.txt && GLIMPSE_concordance_v2.0.0 --gt-val --ac-bins 1 5 10 20 50 100 200 500 1000  2000 5000 10000  20000 50000 100000 146754 --thread 2 --input concordance.txt --output imputed_1k_rp_${METB}_chr${CHR}_rp_af_binning2 && rm -f concordance.txt" --instance-type mem1_ssd1_v2_x2 --priority low --name ${METB}_conc_wgs --folder "${OFOLDER}" -y

	#dx run app-swiss-army-knife -iimage_file="/docker/glimpse2_20221005_73decec-dirty.tar.gz" -icmd="echo \"chr20 /mnt/project/Phasing/PhasingWGS/step8_imputation/reference_panels/sites/rp_s1k_shapeit5_chr20.bcf /mnt/project/data/ukb_wgs/phased_shapeit/panels/validation_greml/wgs_ukb23352_c20_qc_v1_all.bcf /mnt/project/Phasing/PhasingWGS/step8_imputation/ligated/imputed_1k_rp_${METS}_chr${CHR}.bcf\" > concordance.txt && GLIMPSE_concordance_v2.0.0 --gt-val --ac-bins 1 5 10 20 50 100 200 500 1000  2000 5000 10000  20000 50000 100000 146754 --thread 2 --input concordance.txt --output imputed_1k_rp_${METS}_chr${CHR}_rp_af_binning2 && rm -f concordance.txt" --instance-type mem1_ssd1_v2_x2 --priority low --name ${METS}_conc_wgs --folder "${OFOLDER}" -y

	#dx run app-swiss-army-knife -iimage_file="/docker/glimpse2_20221005_73decec-dirty.tar.gz" -icmd="echo \"chr20 /mnt/project/Phasing/PhasingWGS/step8_imputation/reference_panels/sites/rp_s1k_shapeit5_chr20.bcf /mnt/project/data/ukb_wgs/phased_shapeit/panels/validation_greml/wgs_ukb23352_c20_qc_v1_all.bcf /mnt/project/Phasing/PhasingWGS/step8_imputation/ligated/imputed_1k_rp_hrc_chr${CHR}.bcf\" > concordance.txt && GLIMPSE_concordance_v2.0.0 --gt-val --ac-bins 1 5 10 20 50 100 200 500 1000  2000 5000 10000  20000 50000 100000 146754 --thread 2 --input concordance.txt --output imputed_1k_rp_hrc_chr${CHR}_rp_af_binning2 && rm -f concordance.txt" --instance-type mem1_ssd1_v2_x2 --priority low --name hrc_conc_wgs --folder "${OFOLDER}" -y

	
	JOBID2=$(dx run app-swiss-army-knife -iimage_file="/docker/glimpse_v2.0.0-27-g0919952_20221207.tar.gz" -icmd="echo \"chr20 /mnt/project/Phasing/PhasingWGS/step8_imputation/reference_panels/sites/chr20/rp_ukb23352_c20_qc_v1_all_sites.bcf /mnt/project/data/ukb_wgs/phased_shapeit/panels/validation_greml/wgs_ukb23352_c20_qc_v1_all.bcf /mnt/project/PhasingWGS_Official_release/step2_phase_rare/imputation/imputed_1k_rp_shapeit5_chr20.vcf.gz\" > concordance.txt && GLIMPSE2_concordance --gt-val --ac-bins 1 5 10 20 50 100 200 500 1000 2000 5000 10000  20000 50000 100000 149754 --thread 4 --input concordance.txt --output imputed_1k_rp_shapeit5_chr20 && rm -f concordance.txt" --instance-type mem2_ssd1_v2_x4 --priority low --name cnc_rp_shapeit5 --folder "/PhasingWGS_Official_release/step2_phase_rare/concordance/" -y | tail -n1 | cut -d" " -f3)


	JOBID2=$(dx run app-swiss-army-knife -iimage_file="/docker/glimpse_v2.0.0-27-g0919952_20221207.tar.gz" -icmd="echo \"chr20 /mnt/project/Phasing/PhasingWGS/step8_imputation/reference_panels/sites/chr20/rp_ukb23352_c20_qc_v1_all_sites.bcf /mnt/project/data/ukb_wgs/phased_shapeit/panels/validation_greml/wgs_ukb23352_c20_qc_v1_all.bcf /mnt/project/PhasingWGS_Official_release/step3_phase_beagle/imputation/imputed_1k_rp_beagle5.4_chr20.vcf.gz\" > concordance.txt && GLIMPSE2_concordance --gt-val --ac-bins 1 5 10 20 50 100 200 500 1000 2000 5000 10000  20000 50000 100000 149754 --thread 4 --input concordance.txt --output imputed_1k_rp_beagle5.4_chr20 && rm -f concordance.txt" --instance-type mem2_ssd1_v2_x4 --priority low --name cnc_rp_beagle5.4 --folder "/PhasingWGS_Official_release/step3_phase_beagle/concordance/" -y | tail -n1 | cut -d" " -f3)

