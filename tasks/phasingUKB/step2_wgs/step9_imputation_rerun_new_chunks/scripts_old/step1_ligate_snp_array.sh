#!/bin/bash

CHR=chr20
OFOLDER=/Phasing/PhasingWGS/step8_imputation/imputation
METB="beagle5.3"
METS="shapeit5"

#for MET in beagle5.3 shapeit5; do
for MET in shapeit5; do
	OFOLDER=/Phasing/PhasingWGS/step8_imputation/imputation	
	#dx run app-swiss-army-knife -iimage_file="docker/shapeit_ligate_0.0.1.tar.gz" -icmd="ls -1v /mnt/project/Phasing/PhasingWGS/step8_imputation/imputation/${MET}/imputed_100_rp_${MET}_*.vcf.gz > list.txt && shapeit_ligate_v0.0.1 --input list.txt --output imputed_100_rp_${MET}_${CHR}.vcf.gz --index --thread 2 && bcftools index -f imputed_100_rp_${MET}_${CHR}.vcf.gz --threads 2 && rm list.txt" --tag ligate --instance-type mem2_ssd1_v2_x2 --priority normal --name ${MET}_ligate --folder "${OFOLDER}" -y

	OFOLDER=/UKB_PHASING_WGS/step8_imputation/imputation
	dx run app-swiss-army-knife -iimage_file="docker/shapeit_ligate_0.0.1.tar.gz" -icmd="ls -1v /mnt/project/UKB_PHASING_WGS/step8_imputation/imputation/${MET}/imputed_100_rp_${MET}_*.vcf.gz > list.txt && shapeit_ligate_v0.0.1 --input list.txt --output imputed_100_rp_${MET}_${CHR}.vcf.gz --index --thread 2 && bcftools index -f imputed_100_rp_${MET}_${CHR}.vcf.gz --threads 2 && rm list.txt" --tag ligate --instance-type mem2_ssd1_v2_x2 --priority normal --name ${MET}_ligate --folder "${OFOLDER}" -y
done

