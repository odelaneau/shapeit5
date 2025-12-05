#!/bin/bash

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
		printf -v IDG "%03d" ${CNK}
	
		ONAME=HRC.hg38.chr${CHR}
	
		JOBID2=$(dx run app-swiss-army-knife -iin="/docker/beagle.19Apr22.7c0.jar" -icmd="java -Xmx16G -jar beagle.19Apr22.7c0.jar gt=/mnt/project/Phasing/PhasingWGS/step8_imputation/target_snp_array/axiom_1k_c${CHR}_b0_v2.b38.vcf.gz map=/mnt/project/data/plink_maps/plink.prefix.chr20.GRCh38.map ref=/mnt/project/data/HRC/bref/${ONAME}_chr${CHR}_${REGS}_${REGE}.bref3 out=imputed_1k_rp_hrc_chr${CHR}_${REGS}_${REGE} window=500.0  nthreads=4 chrom=${IRG} impute=true gp=true && bcftools index -f imputed_1k_rp_hrc_chr${CHR}_${REGS}_${REGE}.vcf.gz --threads 4" --folder="/Phasing/PhasingWGS/step8_imputation/imputation/hrc" --instance-type mem2_ssd1_v2_x4 --priority low --name imp_rp_hrc_${IDG} -y | tail -n1 | cut -d" " -f3)

		NLINE=$((NLINE+1))
	done < data/new_chr${CHR}_chunks.txt
done
