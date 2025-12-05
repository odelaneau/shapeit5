#!bin/bash

BIN=docker/shapeit5_v1.0.0.tar.gz
ODIR=PhasingWGS_Official_release/step1_phase_common/

threads=36
for CHR in {1..22}; do
	dx run app-swiss-army-knife -iimage_file="/${BIN}" --folder="${ODIR}/" -icmd="ls -1v /mnt/project/PhasingWGS_Official_release/step1_phase_common/chunks/UKB_chr${CHR}.chunk_*.shapeit5_common.bcf > list_ligate.chr${CHR}.txt && SHAPEIT5_ligate_static_v1.0.0 --input list_ligate.chr${CHR}.txt --output UKB_chr${CHR}.shapeit5_common_ligate.bcf --thread ${threads} --index" --instance-type mem1_ssd1_v2_x8 --priority low --name WGS_shapeit5_chr${CHR}_ligate -y

done

