#!bin/bash

for CHR in X; do

	ODIR0=QC/chr${CHR}/support
	dx mkdir -p ${ODIR0}


	SDIR=QC/support
	# 1. upload script
	SCR=script_stats_qc.chrX.sh
	dx upload ${SCR} --path="${SDIR}/"
	
	# 2. upload sample file
	SAMP=ukb_200k.samples_HWE_females.txt # similar as for the autosomes, we computed population metrics on caucasians only. Specifically, on chromosome X, we used only females. This file contains two columns: col1=sample id; col2=population tag. For col2, we used CAU_F for caucasian females, CAU_M for caucasian males, UN_F for the non-caucasian females and UN_M for the non-caucasian males.

	N=$(($(ls ../autosomes/chunks_qc/Splitted/chr${CHR}/chr${CHR}_chunks* | wc -l)-1))



	for SPL in $(seq 0 1 3); do
		echo ${CHR}; echo ${SPL}


		printf -v ID "%03d" $(echo $SPL)
	
		# 3. upload chunk file
		PARAM=../chunks_qc/Splitted/chr${CHR}/chr${CHR}_chunks${ID}
		dx rm "${ODIR0}/$(basename ${PARAM})"
		dx upload ${PARAM} --path="${ODIR0}/"
		PARAM=$(basename ${PARAM})
	
		# 4. run script in parallele
		ODIR1=QC/chr${CHR}/chunks
		dx mkdir -p ${ODIR1}
	
		# 4.1 concat the chunks.
		IN1=ukb24304_chr${CHR}_b*_v1.qceed.bcf
		OUT1=ukb24304_chr${CHR}_chunk_${SPL}.qceed.bcf
		IN2=ukb24304_chr${CHR}_b*_v1.NOqceed.stats.bcf
		OUT2=ukb24304_chr${CHR}_chunk_${SPL}.NOqceed.stats.bcf
		threads=16
		
		dx run app-swiss-army-knife -iin="${ODIR0}/${PARAM}" -iin="${SDIR}/${SCR}" -iin="${SDIR}/${SAMP}" -icmd="cat ${PARAM} | xargs -P 8 -n2 bash ${SCR} && ls -1v ${IN1} > concat_list_chr${CHR}.txt && bcftools concat --threads ${threads} --naive-force -f concat_list_chr${CHR}.txt -o ${OUT1} && bcftools index ${OUT1} --threads ${threads} && rm concat_list_chr${CHR}.txt && rm ${IN1}* && ls -1v ${IN2} > concat_list2_chr${CHR}.txt && bcftools concat --threads ${threads} --naive-force -f concat_list2_chr${CHR}.txt -o ${OUT2} && bcftools index ${OUT2} --threads ${threads} && rm concat_list2_chr${CHR}.txt && rm ${IN2}*" --instance-type mem2_ssd1_v2_x16 --folder="${ODIR1}" --name QC_chr${CHR}_spl_${SPL} --priority normal -y --brief --ignore-reuse

	
	done
done












