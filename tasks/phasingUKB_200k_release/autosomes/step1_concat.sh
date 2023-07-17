#!bin/bash


# This script concat the QC outputs into chromosome-wide files.
#We observed that when using a large number of files as input, the mounted folder that can be accesses using "/mnt/project/" has some problem. We therefore use the -iin command to specifiy the input files. For a large number of input files, we built a script to write the command.


mkdir -p TMP

CHR_START=$1
CHR_END=$2



for CHR in $(seq $CHR_START 1 $CHR_END); do

	# QC concat
	SCR=TMP/script_concat_chr${CHR}.sh
	ODIR=QC/chr${CHR}
	mkdir -p ${ODIR}
	OUT=ukb24304_chr${CHR}.qceed.bcf
	threads=36
	N=$(($(ls chunks_qc/Splitted/chr${CHR}/chr${CHR}_chunks* | wc -l)-1))
	printf '#!bin/bash\ndx run app-swiss-army-knife ' > ${SCR}
	LINE=''
	IDIR=QC/chr${CHR}/chunks
	for SPL in $(seq 0 1 ${N}); do
		LINE="$LINE"-iin="${IDIR}"/ukb24304_chr"${CHR}"_chunk_"${SPL}".qceed.bcf" "		
	done
	LINE2='-icmd="ls -1v ukb24304_chr'${CHR}'_chunk_*.qceed.bcf > concat_list_chr'${CHR}'.txt && bcftools concat --threads '${threads}' --naive-force -f concat_list_chr'${CHR}'.txt -o '${OUT}' && bcftools index '${OUT}' --threads '${threads}' && rm concat_list_chr'${CHR}'.txt" --instance-type mem1_ssd1_v2_x36 --folder='${ODIR}' --name concat_qc_chr'${CHR}' --priority normal -y'
	L="$LINE""$LINE2"
	echo ${L} >> ${SCR}
	bash ${SCR}


	# Stats concat
	SCR=TMP/script_concat_chr${CHR}.stats.sh
	ODIR=QC/Stats/chr${CHR}
	mkdir -p ${ODIR}
	OUT=ukb24304_chr${CHR}.NOqc.stats.bcf
	threads=36
	N=$(($(ls chunks_qc/Splitted/chr${CHR}/chr${CHR}_chunks* | wc -l)-1))
	printf '#!bin/bash\ndx run app-swiss-army-knife ' > ${SCR}
	LINE=''
	IDIR=QC/chr${CHR}/chunks
	for SPL in $(seq 0 1 ${N}); do
		LINE="$LINE"-iin="${IDIR}"/ukb24304_chr"${CHR}"_chunk_"${SPL}".NOqceed.stats.bcf" "		
	done
	LINE2='-icmd="ls -1v ukb24304_chr'${CHR}'_chunk_*.NOqceed.stats.bcf > concat_list_chr'${CHR}'.txt && bcftools concat --threads '${threads}' --naive-force -f concat_list_chr'${CHR}'.txt -o '${OUT}' && bcftools index '${OUT}' --threads '${threads}' && rm concat_list_chr'${CHR}'.txt" --instance-type mem1_ssd1_v2_x36 --folder='${ODIR}' --name concat_stats_chr'${CHR}' --priority normal -y'
	L="$LINE""$LINE2"
	echo ${L} >> ${SCR}
	bash ${SCR}


done
