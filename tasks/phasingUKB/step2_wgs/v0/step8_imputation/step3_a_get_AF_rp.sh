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

		#Setting 2

		JOBID0=$(dx run app-swiss-army-knife -icmd="bcftools view /mnt/project/Phasing/PhasingWGS/step8_imputation/reference_panels/rp_s1k_${METS}_chr${CHR}_${REGS}_${REGE}.vcf.gz -G -r ${ORG} -t ${ORG} -Ob -o rp_s1k_${METS}_chr${CHR}_${OREGS}_${OREGE}.bcf --threads 4 && bcftools index -f rp_s1k_${METS}_chr${CHR}_${OREGS}_${OREGE}.bcf --threads 4" --folder="/Phasing/PhasingWGS/step8_imputation/reference_panels/sites" --instance-type mem1_ssd1_v2_x4 --priority low --name rp_pos_${IDG} -y | tail -n1 | cut -d" " -f3) 
		NLINE=$((NLINE+1))
	done < ../step2_splitchunks/chr20.size4Mb.txt
done
