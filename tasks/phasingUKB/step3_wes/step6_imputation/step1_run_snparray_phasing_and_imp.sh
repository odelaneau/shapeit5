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

for CHR in 13 18; do
	####################################
	#Setting 1
	ONAME="UKB_chr${CHR}.exome_array.WO_parents"
	#JOBID0=$(dx run app-swiss-army-knife -iin="/Phasing/PhasingWGS/step8_imputation/support/1k_samples_wbi.txt" -icmd="bcftools view /mnt/project/${PATHB}/chr${CHR}/${ONAME}.${METB2}.${EXTB} -S ^1k_samples_wbi.txt -Ou  --threads 4 | bcftools filter -e 'AC==0 || AC==AN' -Oz -o rp_s1k_${METB}_chr${CHR}.vcf.gz --threads 4 && bcftools index -f rp_s1k_${METB}_chr${CHR}.vcf.gz --threads 4" --folder="/UKB_PHASING_EXOME_ARRAY/step6_imputation/reference_panels" --instance-type mem1_ssd1_v2_x4 --priority low --name subs_${METB}_${CHR} -y | tail -n1 | cut -d" " -f3) 
	
	##beagle imp
	#JOBID1b=$(dx run app-swiss-army-knife -iin="/docker/bref3.19Apr22.7c0.jar" -icmd="java -Xmx16G -jar bref3.19Apr22.7c0.jar /mnt/project/UKB_PHASING_EXOME_ARRAY/step6_imputation/reference_panels/rp_s1k_${METB}_chr${CHR}.vcf.gz > rp_s1k_${METB}_chr${CHR}.bref3" --folder="/UKB_PHASING_EXOME_ARRAY/step6_imputation/reference_panels/bref/" --instance-type mem3_ssd1_v2_x2 --priority low --name bref_${CHR}_${METB} -y | tail -n1 | cut -d" " -f3) #--depends-on ${JOBID0}
	#JOBID2=$(dx run app-swiss-army-knife -iin="/docker/beagle.19Apr22.7c0.jar" -icmd="java -Xmx120G -jar beagle.19Apr22.7c0.jar gt=/mnt/project/Phasing/PhasingWGS/step8_imputation/target_snp_array/axiom_1k_c${CHR}_b0_v2.b38.vcf.gz map=/mnt/project/data/plink_maps/plink.prefix.chr${CHR}.GRCh38.map ref=/mnt/project/UKB_PHASING_EXOME_ARRAY/step6_imputation/reference_panels/bref/rp_s1k_${METB}_chr${CHR}.bref3 out=imputed_1k_rp_${METB}_chr${CHR} nthreads=10 impute=true gp=true && bcftools index -f imputed_1k_rp_${METB}_chr${CHR}.vcf.gz --threads 4" --folder="/UKB_PHASING_EXOME_ARRAY/step6_imputation/imputation/${METB}" --instance-type mem3_ssd1_v2_x16 --priority low --name imp_${METB}_${CHR} -y | tail -n1 | cut -d" " -f3) #--depends-on ${JOBID0} 

	#from vcf.gz		
	#JOBID2=$(dx run app-swiss-army-knife -iin="/docker/beagle.19Apr22.7c0.jar" -icmd="java -Xmx120G -jar beagle.19Apr22.7c0.jar gt=/mnt/project/Phasing/PhasingWGS/step8_imputation/target_snp_array/axiom_1k_c${CHR}_b0_v2.b38.vcf.gz map=/mnt/project/data/plink_maps/plink.prefix.chr${CHR}.GRCh38.map ref=/mnt/project/UKB_PHASING_EXOME_ARRAY/step6_imputation/reference_panels/rp_s1k_${METB}_chr${CHR}.vcf.gz out=imputed_1k_rp_${METB}_chr${CHR} nthreads=8 impute=true gp=true && bcftools index -f imputed_1k_rp_${METB}_chr${CHR}.vcf.gz --threads 4" --folder="/UKB_PHASING_EXOME_ARRAY/step6_imputation/imputation/${METB}" --instance-type mem3_ssd1_v2_x16 --priority low --name imp_${METB}_${CHR} -y | tail -n1 | cut -d" " -f3) #--depends-on ${JOBID0} 

	####################################
	#Setting 2
	ONAME="UKB_chr${CHR}.exome_array.full"
	#JOBID0=$(dx run app-swiss-army-knife -iin="/Phasing/PhasingWGS/step8_imputation/support/1k_samples_wbi.txt" -icmd="bcftools view /mnt/project/${PATHS}/${ONAME}.${METS}.${EXTS} -S ^1k_samples_wbi.txt -Ou  --threads 4 | bcftools filter -e 'AC==0 || AC==AN' -Oz -o rp_s1k_${METS}_chr${CHR}.vcf.gz --threads 4 && bcftools index -f rp_s1k_${METS}_chr${CHR}.vcf.gz --threads 4" --folder="/UKB_PHASING_EXOME_ARRAY/step6_imputation/reference_panels" --instance-type mem1_ssd1_v2_x4 --priority low --name subs_${METS}_${CHR} -y | tail -n1 | cut -d" " -f3) 

	##shapeit imp
	#JOBID1b=$(dx run app-swiss-army-knife -iin="/docker/bref3.19Apr22.7c0.jar" -icmd="java -Xmx32G -jar bref3.19Apr22.7c0.jar /mnt/project/UKB_PHASING_EXOME_ARRAY/step6_imputation/reference_panels/rp_s1k_${METS}_chr${CHR}.vcf.gz > rp_s1k_${METS}_chr${CHR}.bref3" --folder="/UKB_PHASING_EXOME_ARRAY/step6_imputation/reference_panels/bref/" --instance-type mem3_ssd1_v2_x4 --priority low --name bref_${CHR}_${METS} -y | tail -n1 | cut -d" " -f3) #mem3_ssd1_v2_x2
	JOBID2=$(dx run app-swiss-army-knife -iin="/docker/beagle.19Apr22.7c0.jar" -icmd="java -Xmx120G -jar beagle.19Apr22.7c0.jar gt=/mnt/project/Phasing/PhasingWGS/step8_imputation/target_snp_array/axiom_1k_c${CHR}_b0_v2.b38.vcf.gz map=/mnt/project/data/plink_maps/plink.prefix.chr${CHR}.GRCh38.map ref=/mnt/project/UKB_PHASING_EXOME_ARRAY/step6_imputation/reference_panels/bref/rp_s1k_${METS}_chr${CHR}.bref3 out=imputed_1k_rp_${METS}_chr${CHR} nthreads=10 impute=true gp=true && bcftools index -f imputed_1k_rp_${METS}_chr${CHR}.vcf.gz --threads 4" --folder="/UKB_PHASING_EXOME_ARRAY/step6_imputation/imputation/${METS}" --instance-type mem3_ssd1_v2_x16 --priority low --name imp_${METS}_${CHR} -y | tail -n1 | cut -d" " -f3) #--depends-on ${JOBID0}
	#from vcf.gz		
	#JOBID2=$(dx run app-swiss-army-knife -iin="/docker/beagle.19Apr22.7c0.jar" -icmd="java -Xmx120G -jar beagle.19Apr22.7c0.jar gt=/mnt/project/Phasing/PhasingWGS/step8_imputation/target_snp_array/axiom_1k_c${CHR}_b0_v2.b38.vcf.gz map=/mnt/project/data/plink_maps/plink.prefix.chr${CHR}.GRCh38.map ref=/mnt/project/UKB_PHASING_EXOME_ARRAY/step6_imputation/reference_panels/rp_s1k_${METS}_chr${CHR}.vcf.gz out=imputed_1k_rp_${METS}_chr${CHR} nthreads=8 impute=true gp=true && bcftools index -f imputed_1k_rp_${METS}_chr${CHR}.vcf.gz --threads 4" --folder="/UKB_PHASING_EXOME_ARRAY/step6_imputation/imputation/${METS}" --instance-type mem3_ssd1_v2_x16 --priority low --name imp_${METS}_${CHR} -y | tail -n1 | cut -d" " -f3) #--depends-on ${JOBID0}

done
