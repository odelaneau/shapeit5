#!bin/bash


threads=8


for CHR in {1..22}; do

	ODIR=UKB_PHASING_EXOME_ARRAY/step1_chunks/chr${CHR}
	dx mkdir -p ${ODIR}
	
	
	
	# step1. get sites
	IN=/mnt/project/UKB_PHASING_EXOME_ARRAY/step0_merge/chr${CHR}/UKB.chr${CHR}.exome_array.full.sorted.bcf
	OUT=UKB.chr${CHR}.exome_array.sorted.sites.bcf
	dx run app-swiss-army-knife -icmd="bcftools view --threads ${threads} -G -Ob -o  ${OUT} ${IN} && bcftools index -f ${OUT} --threads ${threads}" --tag sites --instance-type mem1_ssd1_v2_x8 --folder="./${ODIR}/" --name get_sites --priority normal -y


	
	# step2. chunks
	IN=~/dxfuse/DIR/PhasingUKB/UKB_PHASING_EXOME_ARRAY/step1_chunks/chr${CHR}/UKB.chr${CHR}.exome_array.sorted.sites.bcf
	OUT=chunks.chr${CHR}.txt
	#./GLIMPSE_v1.1.1_chunk_static --input ${IN} --window-size 2000000 --buffer-size 250000 --region chr${CHR} --output ${OUT}
	#./GLIMPSE_v1.1.1_chunk_static --input ${IN} --window-size 800000 --window-count 30000 --buffer-size 250000 --region chr${CHR} --output ${OUT}

done



