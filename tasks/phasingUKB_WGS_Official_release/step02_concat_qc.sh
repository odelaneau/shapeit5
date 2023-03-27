#!bin/bash

for CHR in 20; do
	ODIR1=PhasingWGS_Official_release/step0_qc/chunks
	ODIR2=PhasingWGS_Official_release/step0_qc/support
	ODIR3=PhasingWGS_Official_release/step0_qc
	dx mkdir -p ${ODIR2}
	dx mkdir -p ${ODIR3}
	
	dx find data --folder "${ODIR1}" --name "*c${CHR}*.bcf" --delim | sort -k 4 -t$'\t' -V | awk '{print "/mnt/project"$6}' > concat_chr${CHR}.txt
	dx upload concat_chr${CHR}.txt --path="${ODIR2}";
	JOBID0=$(dx run app-swiss-army-knife -icmd="bcftools concat -f /mnt/project/${ODIR2}/concat_chr${CHR}.txt -n -o ukb23352_c${CHR}_qc_v1.bcf --threads 2 && bcftools index ukb23352_c${CHR}_qc_v1.bcf --threads 2 -f" --instance-type mem1_ssd1_v2_x2 --folder="${ODIR3}" --name concat_qc --priority normal -y | tail -n1 | cut -d" " -f3)
	JOBID1=$(dx run app-swiss-army-knife -icmd="bcftools view -G /mnt/project/${ODIR3}/ukb23352_c${CHR}_qc_v1.bcf --threads 4 -Ob -o ukb23352_c${CHR}_qc_v1_sites.bcf && bcftools index ukb23352_c${CHR}_qc_v1_sites.bcf -f --threads 4" --tag sites --instance-type mem1_ssd1_v2_x4 --folder="${ODIR3}" --name get_sites --priority normal --depends-on ${JOBID0} -y | tail -n1 | cut -d" " -f3)
done
