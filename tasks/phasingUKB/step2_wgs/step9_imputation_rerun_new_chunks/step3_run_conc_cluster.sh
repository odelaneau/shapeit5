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

NLINE=0
CNK=0
while read LINE; do 
	CNK=$(echo ${LINE} | cut -d" " -f1)
	CHR2=$(echo ${LINE} | cut -d" " -f2)
	IRG=$(echo ${LINE} | cut -d" " -f3)
	ORG=$(echo ${LINE} | cut -d" " -f4)
	REGS=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f1)
	REGE=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f2)
	OREGS=$(echo ${ORG} | cut -d":" -f 2 | cut -d"-" -f1)
	OREGE=$(echo ${ORG} | cut -d":" -f 2 | cut -d"-" -f2)

	printf -v IDG "%03d" ${CNK}
	
	CHR=$(echo "${CHR2}" | sed -e "s/^chr//")
	####################################

	#Setting 2PhasingWGS_Official_release/step2_phase_rare/concat_benchmark/
	#JOBID0=$(dx run app-swiss-army-knife -iin="/Phasing/PhasingWGS/step8_imputation/support/1k_samples_wbi.txt" -icmd="bcftools view /mnt/project/PhasingWGS_Official_release/step2_phase_rare/concat_benchmark/${CHR2}/UKB_${CHR2}.full.shapeit5_rare.bcf -r ${IRG} -S ^1k_samples_wbi.txt -Ou  --threads 4 | bcftools filter -e 'AC==0 || AC==AN' -Oz -o UKB_${CHR2}_${REGS}_${REGE}.full.shapeit5_rare.vcf.gz --threads 4 && bcftools index -f UKB_${CHR2}_${REGS}_${REGE}.full.shapeit5_rare.vcf.gz --threads 4" --folder="/PhasingWGS_Official_release/step2_phase_rare/concat_benchmark/${CHR2}/reference_panels" --instance-type mem1_ssd1_v2_x4 --priority low --name subset_${METS}_${IDG} -y | tail -n1 | cut -d" " -f3) 

	##shapeit imp
	#from vcf.gz		
	JOBID2=$(dx run app-swiss-army-knife -iimage_file="/docker/glimpse_v2.0.0-27-g0919952_20221207.tar.gz" -icmd="echo \"${ORG} /mnt/project/Phasing/PhasingWGS/step8_imputation/reference_panels/sites/${CHR2}/rp_ukb23352_c${CHR}_qc_v1_all_sites.bcf /mnt/project/data/ukb_wgs/phased_shapeit/panels/validation_greml/wgs_ukb23352_c${CHR}_qc_v1_all.bcf /mnt/project/PhasingWGS_Official_release/step2_phase_rare/concat_benchmark/${CHR2}/imputation/imputed_1k_rp_shapeit5_${CHR2}_${REGS}_${REGE}.vcf.gz\" > concordance.txt && GLIMPSE2_concordance --gt-val --ac-bins 1 5 10 20 50 100 200 500 1000  2000 5000 10000  20000 50000 100000 149754 --thread 2 --input concordance.txt --output imputed_1k_rp_shapeit5_${CHR2}_${REGS}_${REGE} && rm -f concordance.txt" --instance-type mem1_ssd1_v2_x2 --priority low --name ${CHR} --folder "/PhasingWGS_Official_release/step2_phase_rare/concat_benchmark/${CHR2}/concordance/" -y | tail -n1 | cut -d" " -f3)

	NLINE=$((NLINE+1))
done < data/new.chunks.rerun.size4Mb.txt
