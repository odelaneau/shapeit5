#!/bin/bash

REF=HRC
#/mnt/project/data/ukb_wgs/phased_beagle/panels/reference/sites/ukb_white_british_sites.bcf
for COV in 0.1 1.0 4.0; do
	CHR=20
	MET=glimpse2
	OFOLDER=/test/${REF}/concordance_test
	#dx run app-swiss-army-knife -iimage_file="/docker/glimpse2_0.0.1.tar.gz" -icmd="echo \"chr20 /mnt/project/data/HRC/HRC.hg38.chr20_sites.bcf /mnt/project/data/ukb_wgs/phased_beagle/panels/validation/val_ukb23352_c20_qc_v1_all.bcf /mnt/project/test/${REF}/imputed_c${CHR}_${MET}_rp${REF}_${COV}x_100.bcf\" > concordance.txt && GLIMPSE_concordance_v2.0.0 --gt-val --allele-counts --thread 2 --input concordance.txt --output c${CHR}_${MET}_rp${REF}_${COV}x_100 && rm -f concordance.txt" --tag concordance --instance-type mem1_ssd1_v2_x2 --priority normal --tag cov${COV}x --name ${MET}_concordance --folder "${OFOLDER}" -y
	#dx run app-swiss-army-knife -iimage_file="/docker/glimpse2_0.0.1.tar.gz" -icmd="echo \"chr20 /mnt/project/data/HRC/HRC.hg38.chr20_sites_isec_quilt.bcf /mnt/project/data/ukb_wgs/phased_beagle/panels/validation/val_ukb23352_c20_qc_v1_all.bcf /mnt/project/test/${REF}/imputed_c${CHR}_${MET}_rp${REF}_${COV}x_100.bcf\" > concordance.txt && GLIMPSE_concordance_v2.0.0 --gt-val --allele-counts --thread 2 --input concordance.txt --output c${CHR}_${MET}_rp${REF}_${COV}x_100_isec && rm -f concordance.txt" --tag concordance --instance-type mem1_ssd1_v2_x2 --priority normal --tag cov${COV}x --name ${MET}_concordance --folder "${OFOLDER}" -y

	MET=glimpse1.1.1
	#dx run app-swiss-army-knife -iimage_file="/docker/glimpse2_0.0.1.tar.gz" -icmd="echo \"chr20 /mnt/project/data/HRC/HRC.hg38.chr20_sites.bcf /mnt/project/data/ukb_wgs/phased_beagle/panels/validation/val_ukb23352_c20_qc_v1_all.bcf /mnt/project/test/${REF}/imputed_c${CHR}_${MET}_rp${REF}_${COV}x_100.bcf\" > concordance.txt && GLIMPSE_concordance_v2.0.0 --gt-val --allele-counts --thread 2 --input concordance.txt --output c${CHR}_${MET}_rp${REF}_${COV}x_100 && rm -f concordance.txt" --tag concordance --instance-type mem1_ssd1_v2_x2 --priority normal --tag cov${COV}x --name ${MET}_concordance --folder "${OFOLDER}" -y
	dx run app-swiss-army-knife -iimage_file="/docker/glimpse2_0.0.1.tar.gz" -icmd="echo \"chr20 /mnt/project/data/HRC/HRC.hg38.chr20_sites_isec_quilt.bcf /mnt/project/data/ukb_wgs/phased_beagle/panels/validation/val_ukb23352_c20_qc_v1_all.bcf /mnt/project/test/${REF}/imputed_c${CHR}_${MET}_rp${REF}_${COV}x_100.bcf\" > concordance.txt && GLIMPSE_concordance_v2.0.0 --gt-val --allele-counts --thread 2 --input concordance.txt --output c${CHR}_${MET}_rp${REF}_${COV}x_100_isec && rm -f concordance.txt" --tag concordance --instance-type mem1_ssd1_v2_x2 --priority normal --tag cov${COV}x --name ${MET}_concordance --folder "${OFOLDER}" -y

	MET=quilt
	#dx run app-swiss-army-knife -iimage_file="/docker/glimpse2_0.0.1.tar.gz" -icmd="echo \"chr20 /mnt/project/data/HRC/HRC.hg38.chr20_sites.bcf /mnt/project/data/ukb_wgs/phased_beagle/panels/validation/val_ukb23352_c20_qc_v1_all.bcf /mnt/project/test/${REF}/imputed_c${CHR}_${MET}_rp${REF}_${COV}x_100.vcf.gz\" > concordance.txt && GLIMPSE_concordance_v2.0.0 --gt-val --allele-counts --thread 2 --input concordance.txt --output c${CHR}_${MET}_rp${REF}_${COV}x_100 && rm -f concordance.txt" --tag concordance --instance-type mem1_ssd1_v2_x2 --priority normal --tag cov${COV}x --name ${MET}_concordance --folder "${OFOLDER}" -y
	#dx run app-swiss-army-knife -iimage_file="/docker/glimpse2_0.0.1.tar.gz" -icmd="echo \"chr20 /mnt/project/data/HRC/HRC.hg38.chr20_sites_isec_quilt.bcf /mnt/project/data/ukb_wgs/phased_beagle/panels/validation/val_ukb23352_c20_qc_v1_all.bcf /mnt/project/test/${REF}/imputed_c${CHR}_${MET}_rp${REF}_${COV}x_100.vcf.gz\" > concordance.txt && GLIMPSE_concordance_v2.0.0 --gt-val --allele-counts --thread 2 --input concordance.txt --output c${CHR}_${MET}_rp${REF}_${COV}x_100_isec && rm -f concordance.txt" --tag concordance --instance-type mem1_ssd1_v2_x2 --priority normal --tag cov${COV}x --name ${MET}_concordance --folder "${OFOLDER}" -y
done

CHR=20
MET=beagle5.4
COV=snp_array
OFOLDER=/test/${REF}/concordance_test
#dx run app-swiss-army-knife -iimage_file="/docker/glimpse2_0.0.1.tar.gz" -icmd="echo \"chr20 /mnt/project/data/HRC/HRC.hg38.chr20_sites.bcf /mnt/project/data/ukb_wgs/phased_beagle/panels/validation/val_ukb23352_c20_qc_v1_all.bcf /mnt/project/test/${REF}/snp_array/imputed_beagle5.4_rp${REF}_snp_array_100_chr20.vcf.gz\" > concordance.txt && GLIMPSE_concordance_v2.0.0 --gt-val --allele-counts --thread 2 --input concordance.txt --output c${CHR}_${MET}_rp${REF}_${COV}x_100 && rm -f concordance.txt" --tag concordance --instance-type mem1_ssd1_v2_x2 --priority normal --tag cov${COV}x --name ${MET}_concordance --folder "${OFOLDER}" -y
#dx run app-swiss-army-knife -iimage_file="/docker/glimpse2_0.0.1.tar.gz" -icmd="echo \"chr20 /mnt/project/data/HRC/HRC.hg38.chr20_sites_isec_quilt.bcf /mnt/project/data/ukb_wgs/phased_beagle/panels/validation/val_ukb23352_c20_qc_v1_all.bcf /mnt/project/test/${REF}/snp_array/imputed_beagle5.4_rp${REF}_snp_array_100_chr20.vcf.gz\" > concordance.txt && GLIMPSE_concordance_v2.0.0 --gt-val --allele-counts --thread 2 --input concordance.txt --output c${CHR}_${MET}_rp${REF}_${COV}x_100_isec && rm -f concordance.txt" --tag concordance --instance-type mem1_ssd1_v2_x2 --priority normal --tag cov${COV}x --name ${MET}_concordance --folder "${OFOLDER}" -y

