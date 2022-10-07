#!/bin/bash

#Setting 1
PATHB="Phasing/PhasingWGS/step3_runbeagle/N147754"
EXTB="vcf.gz"
METB="beagle5.4"

#Setting 2
PATHS="Phasing/PhasingWGS/step5_runshapeit_phase2/N147754"
EXTS="default.bcf"
METS="shapeit5"

for CHR in 20; do
	ONAME="benchmark_ukb23352_c${CHR}_qc_v1"
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

		####################################

		#Setting 1
		
		#JOBID0=$(dx run app-swiss-army-knife -iin="/Phasing/PhasingWGS/step8_imputation/support/1k_samples_wbi.txt" -icmd="bcftools view /mnt/project/${PATHB}/${ONAME}.subset.N147754.${IRG}.${EXTB} -S ^1k_samples_wbi.txt -Ou  --threads 4 | bcftools filter -e 'AC==0 || AC==AN' -Oz -o rp_s1k_${METB}_chr${CHR}_${REGS}_${REGE}.vcf.gz --threads 4 && bcftools index -f rp_s1k_${METB}_chr${CHR}_${REGS}_${REGE}.vcf.gz --threads 4" --folder="/Phasing/PhasingWGS/step8_imputation/reference_panels" --instance-type mem1_ssd1_v2_x4 --priority low --name subset_${METB}_${IDG} -y | tail -n1 | cut -d" " -f3) 
		
		##beagle imp
		#from vcf.gz		
		JOBID2=$(dx run app-swiss-army-knife -iin="/docker/beagle.19Apr22.7c0.jar" -icmd="java -Xmx16G -jar beagle.19Apr22.7c0.jar gt=/mnt/project/Phasing/PhasingWGS/step8_imputation/target_snp_array/axiom_1k_c${CHR}_b0_v2.b38.vcf.gz map=/mnt/project/data/plink_maps/plink.prefix.chr20.GRCh38.map ref=/mnt/project/Phasing/PhasingWGS/step8_imputation/reference_panels/rp_s1k_${METB}_chr${CHR}_${REGS}_${REGE}.vcf.gz out=imputed_1k_rp_${METB}_chr${CHR}_${REGS}_${REGE} window=500.0  nthreads=4 chrom=${IRG} impute=true gp=true && bcftools index -f imputed_1k_rp_${METB}_chr${CHR}_${REGS}_${REGE}.vcf.gz --threads 4" --folder="/Phasing/PhasingWGS/step8_imputation/imputation/${METB}" --instance-type mem3_ssd1_v2_x4 --priority low --name imp_rp_${METB}_${IDG} -y | tail -n1 | cut -d" " -f3) #--depends-on ${JOBID0}

		####################################

		#Setting 2

		#JOBID0=$(dx run app-swiss-army-knife -iin="/Phasing/PhasingWGS/step8_imputation/support/1k_samples_wbi.txt" -icmd="bcftools view /mnt/project/${PATHS}/${ONAME}.subset.N147754.${IRG}.${METS}.${EXTS} -S ^1k_samples_wbi.txt -Ou  --threads 4 | bcftools filter -e 'AC==0 || AC==AN' -Oz -o rp_s1k_${METS}_chr${CHR}_${REGS}_${REGE}.vcf.gz --threads 4 && bcftools index -f rp_s1k_${METS}_chr${CHR}_${REGS}_${REGE}.vcf.gz --threads 4" --folder="/Phasing/PhasingWGS/step8_imputation/reference_panels" --instance-type mem1_ssd1_v2_x4 --priority low --name subset_${METS}_${IDG} -y | tail -n1 | cut -d" " -f3) 

		##shapeit imp
		#from vcf.gz		
		#JOBID2=$(dx run app-swiss-army-knife -iin="/docker/beagle.19Apr22.7c0.jar" -icmd="java -Xmx16G -jar beagle.19Apr22.7c0.jar gt=/mnt/project/Phasing/PhasingWGS/step8_imputation/target_snp_array/axiom_1k_c${CHR}_b0_v2.b38.vcf.gz map=/mnt/project/data/plink_maps/plink.prefix.chr20.GRCh38.map ref=/mnt/project/Phasing/PhasingWGS/step8_imputation/reference_panels/rp_s1k_${METS}_chr${CHR}_${REGS}_${REGE}.vcf.gz out=imputed_1k_rp_${METS}_chr${CHR}_${REGS}_${REGE} window=500.0  nthreads=4 chrom=${IRG} impute=true gp=true && bcftools index -f imputed_1k_rp_${METS}_chr${CHR}_${REGS}_${REGE}.vcf.gz --threads 4" --folder="/Phasing/PhasingWGS/step8_imputation/imputation/${METS}" --instance-type mem2_ssd1_v2_x4 --priority low --name imp_rp_${METS}_${IDG} --depends-on ${JOBID0} -y | tail -n1 | cut -d" " -f3)

		NLINE=$((NLINE+1))
	done < ../step2_splitchunks/chr20.size4Mb.txt
done
