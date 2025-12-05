#!/bin/bash

ONAME="benchmark_ukb23352_c20_qc_v1"

#Setting 1
PATHB="Phasing/PhasingWGS/step3_runbeagle/beagle5.3"
EXTB="vcf.gz.vcf.gz"
METB="beagle5.3"

#Setting 2
PATHS="Phasing/PhasingWGS/step5_runshapeit_phase2/old_beforeSubsetting"
EXTS="scaffold.bcf"
METS="shapeit5"

#Setting Robin
ONAMER="ukb23352_c20_qc_v1"
PATHR="UKB_PHASING_WGS/step3_phasing_rares"
EXTR="rares_0.01.bcf"
METR="shapeit5"

for CHR in 20; do
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
		
		#JOBID0=$(dx run app-swiss-army-knife -iin="/Phasing/PhasingWGS/step8_imputation/target_snp_array/ids_100_samples.txt" -icmd="bcftools view /mnt/project/${PATHB}/${ONAME}.${IRG}.${METB}.${EXTB} -S ^ids_100_samples.txt --force-samples -Ou  --threads 4 | bcftools filter -e 'AC==0 || AC==AN' -Oz -o rp_s100_${METB}_chr${CHR}_${REGS}_${REGE}.vcf.gz --threads 4 && bcftools index -f rp_s100_${METB}_chr${CHR}_${REGS}_${REGE}.vcf.gz --threads 4" --tag rp_subset --folder="/Phasing/PhasingWGS/step8_imputation/reference_panels" --instance-type mem1_ssd1_v2_x4 --priority normal --name subset_${METB}_${IDG} -y | tail -n1 | cut -d" " -f3) 
		
		#Unused ##bref
		#Unused #JOBID1=$(dx run app-swiss-army-knife -iin="/docker/bref3.19Apr22.7c0.jar" -icmd="java -Xmx32G -jar bref3.19Apr22.7c0.jar /mnt/project/Phasing/PhasingWGS/step8_imputation/reference_panels/rp_s100_${METB}_chr${CHR}_${REGS}_${REGE}.vcf.gz > rp_s100_${METB}_chr${CHR}_${REGS}_${REGE}.bref3" --tag bref --folder="/Phasing/PhasingWGS/step8_imputation/reference_panels" --instance-type mem3_ssd1_v2_x4 --priority normal --name bref_${METB}_${IDG} --depends-on ${JOBID0} -y | tail -n1 | cut -d" " -f3) 

		##beagle imp
		#from vcf.gz		
		#JOBID2=$(dx run app-swiss-army-knife -iin="/docker/beagle.19Apr22.7c0.jar" -icmd="java -Xmx16G -jar beagle.19Apr22.7c0.jar gt=/mnt/project/Phasing/PhasingWGS/step8_imputation/target_snp_array/100_samples_c20_b0_v2.b38.r1.vcf.gz map=/mnt/project/data/plink_maps/plink.prefix.chr20.GRCh38.map ref=/mnt/project/Phasing/PhasingWGS/step8_imputation/reference_panels/rp_s100_${METB}_chr${CHR}_${REGS}_${REGE}.vcf.gz out=imputed_100_rp_${METB}_chr${CHR}_${REGS}_${REGE} window=500.0  nthreads=4 chrom=${IRG} impute=true gp=true && bcftools index -f imputed_100_rp_${METB}_chr${CHR}_${REGS}_${REGE}.vcf.gz --threads 4" --tag beagle5.4 --tag imputation --folder="/Phasing/PhasingWGS/step8_imputation/imputation/${METB}" --instance-type mem2_ssd1_v2_x4 --priority normal --name imp_rp_${METS}_${IDG} --depends-on ${JOBID0} -y | tail -n1 | cut -d" " -f3)
	
		#Unused #from bref		
		#Unused #JOBID2=$(dx run app-swiss-army-knife -iin="/docker/beagle.19Apr22.7c0.jar" -icmd="java -Xmx16G -jar beagle.19Apr22.7c0.jar gt=/mnt/project/Phasing/PhasingWGS/step8_imputation/target_snp_array/100_samples_c20_b0_v2.b38.r1.vcf.gz map=/mnt/project/data/plink_maps/plink.prefix.chr20.GRCh38.map ref=/mnt/project/Phasing/PhasingWGS/step8_imputation/reference_panels/rp_s100_${METB}_chr${CHR}_${REGS}_${REGE}.bref3 out=imputed_100_rp_${METB}_chr${CHR}_${REGS}_${REGE} window=500.0  nthreads=4 chrom=${IRG} impute=true gp=true && bcftools index -f imputed_100_rp_${METB}_chr${CHR}_${REGS}_${REGE}.vcf.gz --threads 4" --tag beagle5.4 --tag chr20 --tag imputation --folder="/Phasing/PhasingWGS/step8_imputation/imputation/${METB}" --instance-type mem2_ssd1_v2_x4 --priority normal --name imp_rp_${METS}_${IDG} --depends-on ${JOBID1} -y | tail -n1 | cut -d" " -f3)

		####################################

		#Setting 2

		#JOBID0=$(dx run app-swiss-army-knife -iin="/Phasing/PhasingWGS/step8_imputation/target_snp_array/ids_100_samples.txt" -icmd="bcftools view /mnt/project/${PATHS}/${ONAME}.${IRG}.${METS}.${EXTS} -S ^ids_100_samples.txt --force-samples -Ou  --threads 4 | bcftools filter -e 'AC==0 || AC==AN' -Oz -o rp_s100_${METS}_chr${CHR}_${REGS}_${REGE}.vcf.gz --threads 4 && bcftools index -f rp_s100_${METS}_chr${CHR}_${REGS}_${REGE}.vcf.gz --threads 4" --tag rp_subset --folder="/Phasing/PhasingWGS/step8_imputation/reference_panels" --instance-type mem1_ssd1_v2_x4 --priority normal --name subset_${METS}_${IDG} -y | tail -n1 | cut -d" " -f3) 
		
		#Unused ##bref
		#Unused #JOBID1=$(dx run app-swiss-army-knife -iin="/docker/bref3.19Apr22.7c0.jar" -icmd="java -Xmx32G -jar bref3.19Apr22.7c0.jar /mnt/project/Phasing/PhasingWGS/step8_imputation/reference_panels/rp_s100_${METS}_chr${CHR}_${REGS}_${REGE}.vcf.gz > rp_s100_${METS}_chr${CHR}_${REGS}_${REGE}.bref3" --tag bref --folder="/Phasing/PhasingWGS/step8_imputation/reference_panels" --instance-type mem3_ssd1_v2_x4 --priority normal --name bref_${METS}_${IDG} --depends-on ${JOBID0} -y | tail -n1 | cut -d" " -f3) 

		##beagle imp
		#from vcf.gz		
		#JOBID2=$(dx run app-swiss-army-knife -iin="/docker/beagle.19Apr22.7c0.jar" -icmd="java -Xmx16G -jar beagle.19Apr22.7c0.jar gt=/mnt/project/Phasing/PhasingWGS/step8_imputation/target_snp_array/100_samples_c20_b0_v2.b38.r1.vcf.gz map=/mnt/project/data/plink_maps/plink.prefix.chr20.GRCh38.map ref=/mnt/project/Phasing/PhasingWGS/step8_imputation/reference_panels/rp_s100_${METS}_chr${CHR}_${REGS}_${REGE}.vcf.gz out=imputed_100_rp_${METS}_chr${CHR}_${REGS}_${REGE} window=500.0  nthreads=4 chrom=${IRG} impute=true gp=true && bcftools index -f imputed_100_rp_${METS}_chr${CHR}_${REGS}_${REGE}.vcf.gz --threads 4" --tag beagle5.4 --tag imputation --folder="/Phasing/PhasingWGS/step8_imputation/imputation/${METS}" --instance-type mem2_ssd1_v2_x4 --priority normal --name imp_rp_${METS}_${IDG} --depends-on ${JOBID0} -y | tail -n1 | cut -d" " -f3)
	
		#Unused #from bref		
		#Unused #JOBID2=$(dx run app-swiss-army-knife -iin="/docker/beagle.19Apr22.7c0.jar" -icmd="java -Xmx16G -jar beagle.19Apr22.7c0.jar gt=/mnt/project/Phasing/PhasingWGS/step8_imputation/target_snp_array/100_samples_c20_b0_v2.b38.r1.vcf.gz map=/mnt/project/data/plink_maps/plink.prefix.chr20.GRCh38.map ref=/mnt/project/Phasing/PhasingWGS/step8_imputation/reference_panels/rp_s100_${METS}_chr${CHR}_${REGS}_${REGE}.bref3 out=imputed_100_rp_${METS}_chr${CHR}_${REGS}_${REGE} window=500.0  nthreads=4 chrom=${IRG} impute=true gp=true && bcftools index -f imputed_100_rp_${METS}_chr${CHR}_${REGS}_${REGE}.vcf.gz --threads 4" --tag beagle5.4 --tag chr20 --tag imputation --folder="/Phasing/PhasingWGS/step8_imputation/imputation/${METS}" --instance-type mem2_ssd1_v2_x4 --priority normal --name imp_rp_${METS}_${IDG} --depends-on ${JOBID1} -y | tail -n1 | cut -d" " -f3)

		####################################

		#Setting 3

		#JOBID0=$(dx run app-swiss-army-knife -iin="/Phasing/PhasingWGS/step8_imputation/target_snp_array/ids_100_samples.txt" -icmd="bcftools view /mnt/project/${PATHR}/chr${CHR}/${ONAMER}.${CHR}_${OREGS}_${OREGE}.${METR}.${EXTR} -S ^ids_100_samples.txt --force-samples -Ou  --threads 4 | bcftools filter -e 'AC==0 || AC==AN' -Oz -o rp_s100_${METR}_chr${CHR}_${REGS}_${REGE}.vcf.gz --threads 4 && bcftools index -f rp_s100_${METR}_chr${CHR}_${REGS}_${REGE}.vcf.gz --threads 4" --tag rp_subset --folder="/UKB_PHASING_WGS/step8_imputation/reference_panels" --instance-type mem1_ssd1_v2_x4 --priority normal --name subset_${METS}_${IDG} -y | tail -n1 | cut -d" " -f3) 
		
		#Unused ##bref
		#Unused #JOBID1=$(dx run app-swiss-army-knife -iin="/docker/bref3.19Apr22.7c0.jar" -icmd="java -Xmx32G -jar bref3.19Apr22.7c0.jar /mnt/project/Phasing/PhasingWGS/step8_imputation/reference_panels/rp_s100_${METS}_chr${CHR}_${REGS}_${REGE}.vcf.gz > rp_s100_${METS}_chr${CHR}_${REGS}_${REGE}.bref3" --tag bref --folder="/Phasing/PhasingWGS/step8_imputation/reference_panels" --instance-type mem3_ssd1_v2_x4 --priority normal --name bref_${METS}_${IDG} --depends-on ${JOBID0} -y | tail -n1 | cut -d" " -f3) 

		##beagle imp
		#from vcf.gz		
		JOBID2=$(dx run app-swiss-army-knife -iin="/docker/beagle.19Apr22.7c0.jar" -icmd="java -Xmx16G -jar beagle.19Apr22.7c0.jar gt=/mnt/project/Phasing/PhasingWGS/step8_imputation/target_snp_array/100_samples_c20_b0_v2.b38.r1.vcf.gz map=/mnt/project/data/plink_maps/plink.prefix.chr20.GRCh38.map ref=/mnt/project/UKB_PHASING_WGS/step8_imputation/reference_panels/rp_s100_${METR}_chr${CHR}_${REGS}_${REGE}.vcf.gz out=imputed_100_rp_${METR}_chr${CHR}_${REGS}_${REGE} window=500.0  nthreads=4 chrom=${IRG} impute=true gp=true && bcftools index -f imputed_100_rp_${METR}_chr${CHR}_${REGS}_${REGE}.vcf.gz --threads 4" --tag beagle5.4 --tag imputation --folder="/UKB_PHASING_WGS/step8_imputation/imputation/${METR}" --instance-type mem2_ssd1_v2_x4 --priority normal --name imp_rp_${METS}_${IDG} -y | tail -n1 | cut -d" " -f3)
	
		#Unused #from bref		
		#Unused #JOBID2=$(dx run app-swiss-army-knife -iin="/docker/beagle.19Apr22.7c0.jar" -icmd="java -Xmx16G -jar beagle.19Apr22.7c0.jar gt=/mnt/project/Phasing/PhasingWGS/step8_imputation/target_snp_array/100_samples_c20_b0_v2.b38.r1.vcf.gz map=/mnt/project/data/plink_maps/plink.prefix.chr20.GRCh38.map ref=/mnt/project/Phasing/PhasingWGS/step8_imputation/reference_panels/rp_s100_${METS}_chr${CHR}_${REGS}_${REGE}.bref3 out=imputed_100_rp_${METS}_chr${CHR}_${REGS}_${REGE} window=500.0  nthreads=4 chrom=${IRG} impute=true gp=true && bcftools index -f imputed_100_rp_${METS}_chr${CHR}_${REGS}_${REGE}.vcf.gz --threads 4" --tag beagle5.4 --tag chr20 --tag imputation --folder="/Phasing/PhasingWGS/step8_imputation/imputation/${METS}" --instance-type mem2_ssd1_v2_x4 --priority normal --name imp_rp_${METS}_${IDG} --depends-on ${JOBID1} -y | tail -n1 | cut -d" " -f3)


		NLINE=$((NLINE+1))
	done < ../step2_splitchunks/chr20.size4Mb.txt
done
