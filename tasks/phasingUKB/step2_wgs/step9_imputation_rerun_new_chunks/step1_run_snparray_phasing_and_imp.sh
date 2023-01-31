#!/bin/bash

#Setting 1
PATHB="Phasing/PhasingWGS/step3_runbeagle/N147754"
EXTB="vcf.gz"
METB="beagle5.4"

#Setting 2
PATHS="Phasing/PhasingWGS/step5_runshapeit_phase2/N147754"
EXTS="default.bcf"
METS="shapeit5"

#for CHR in 20; do

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
		JOBID0=$(dx run app-swiss-army-knife -iin="/Phasing/PhasingWGS/step8_imputation/support/1k_samples_wbi.txt" -icmd="bcftools view /mnt/project/PhasingWGS_Official_release/step2_phase_rare/concat_benchmark/${CHR2}/UKB_${CHR2}.full.shapeit5_rare.bcf -r ${IRG} -S ^1k_samples_wbi.txt -Ou  --threads 4 | bcftools filter -e 'AC==0 || AC==AN' -Oz -o UKB_${CHR2}_${REGS}_${REGE}.full.shapeit5_rare.vcf.gz --threads 4 && bcftools index -f UKB_${CHR2}_${REGS}_${REGE}.full.shapeit5_rare.vcf.gz --threads 4" --folder="/PhasingWGS_Official_release/step2_phase_rare/concat_benchmark/${CHR2}/reference_panels" --instance-type mem1_ssd1_v2_x4 --priority low --name subset_${METS}_${IDG} -y | tail -n1 | cut -d" " -f3) 

		##shapeit imp
		#from vcf.gz		
		JOBID2=$(dx run app-swiss-army-knife -iin="/docker/beagle.19Apr22.7c0.jar" -icmd="java -Xmx16G -jar beagle.19Apr22.7c0.jar gt=/mnt/project/Phasing/PhasingWGS/step8_imputation/target_snp_array/axiom_1k_c${CHR}_b0_v2.b38.vcf.gz map=/mnt/project/data/plink_maps/plink.prefix.${CHR2}.GRCh38.map ref=/mnt/project/PhasingWGS_Official_release/step2_phase_rare/concat_benchmark/${CHR2}/reference_panels/UKB_${CHR2}_${REGS}_${REGE}.full.shapeit5_rare.vcf.gz out=imputed_1k_rp_${METS}_chr${CHR}_${REGS}_${REGE} window=500.0  nthreads=4 chrom=${IRG} impute=true gp=true && bcftools index -f imputed_1k_rp_${METS}_chr${CHR}_${REGS}_${REGE}.vcf.gz --threads 4" --folder="/PhasingWGS_Official_release/step2_phase_rare/concat_benchmark/${CHR2}/imputation/" --instance-type mem3_ssd1_v2_x4 --priority low --name imp_rp_${METS}_${IDG} -y --depends-on ${JOBID0} | tail -n1 | cut -d" " -f3)

		NLINE=$((NLINE+1))
	done < data/new.chunks.rerun.size4Mb.txt
#done
