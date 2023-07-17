#!bin/bash


CHR_START=$1
CHR_END=$2



for CHR in $(seq $CHR_START 1 $CHR_END); do

	ODIR0=QC/chr${CHR}/support
	dx mkdir -p ${ODIR0}


	SDIR=/QC/support
	# 1. upload script
	SCR=script_stats_qc.sh
	dx rm ${SDIR}/${SCR}
	dx upload ${SCR} --path="${SDIR}/"
	
	# 2. upload sample file
	SAMP=ukb_200k.samples_with_pop_tag.unrelated_only.txt # this is a two column file with col1= sample_id, col2=population. It contains only unrelated individuals. For col2, we used CAU for caucasians et UN for the remaining samples. It allows us to use caucasian metric only for QC (e.g HWE).

	

	N=$(($(ls chunks_qc/Splitted/chr${CHR}/chr${CHR}_chunks* | wc -l)-1)) # the data before QC is stored in more thank 60k small chunks genome-wide. We groupped these in 100 groups of chunks for each chromosome (since we can run ~100 jobs in parrallel on the RAP). We run on QC job per group of chunks. Within each of these jobs, 8 QC jobs run in parallel (using xargs command). This allows us to considerably speed up the genome-wide running time.


	for SPL in $(seq 0 1 ${N}); do
		echo ${CHR}; echo ${SPL}


		printf -v ID "%03d" $(echo $SPL)
	
		# 3. upload chunk file
		PARAM=chunks_qc/Splitted/chr${CHR}/chr${CHR}_chunks${ID}
#		dx rm "${ODIR0}/$(basename ${PARAM})"
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
