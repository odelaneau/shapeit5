#!bin/bash


PREFIX=ukb200k

mkdir -p TMP

BIN=docker/shapeit5.ukb200k.tar.gz
IDIR=Phasing/step1_phase_common/chunks
ODIR=Phasing/step2_ligate
dx mkdir -p ${ODIR}
threads=8

PED=ukb_200k.ped # sample file as in step2_phase_common.sh

for CHR in 7; do

	SCR=TMP/script_ligate_chr${CHR}.sh
	IN=/mnt/project/${IDIR}/UKB_chr${CHR}.chunk_*.shapeit5_common.bcf
	OUT=${PREFIX}_chr${CHR}.shapeit5_common_ligate.bcf


	printf '#!bin/bash\ndx run app-swiss-army-knife --folder='${ODIR}/' -iimage_file='/${BIN}' ' > ${SCR}
	
	N=$(($(cat chunks/phase_common/chunks_chr${CHR}.txt | wc -l)-1))

	LINE=''

	for SPL in $(seq 0 1 ${N}); do

		LINE="$LINE"-iin="${IDIR}"/"${PREFIX}"_chr"${CHR}".chunk_"${SPL}".shapeit5_common.bcf" "
		LINE="$LINE"-iin="${IDIR}"/"${PREFIX}"_chr"${CHR}".chunk_"${SPL}".shapeit5_common.bcf.csi" "
		
	done

	LINE2='-iin=data/'${PED}' -icmd="ls -1v '${PREFIX}'_chr'${CHR}'.chunk_*.shapeit5_common.bcf > list_ligate.chr'${CHR}'.txt && ligate_static --pedigree '${PED}' --input list_ligate.chr'${CHR}'.txt --output '${OUT}' --thread '${threads}' --index" --instance-type mem1_ssd1_v2_x8 --priority low --name ligate_chr'${CHR}' -y --brief --ignore-reuse'

	L="$LINE""$LINE2"
	echo ${L} >> ${SCR}
	bash ${SCR}

done

