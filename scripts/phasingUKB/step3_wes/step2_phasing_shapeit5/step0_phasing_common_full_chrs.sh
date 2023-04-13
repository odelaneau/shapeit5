#!/bin/bash


for CHR in {1..22}; do


	MAP=/mnt/project/data/shapeit_maps/chr${CHR}.b38.gmap.gz
	
	
	CHUNKS=../step1_chunking/chunks.chr${CHR}.txt
	CUT=0.001
	THREADS=64
	ODIR=UKB_PHASING_EXOME_ARRAY/step2_phasing_common/chr${CHR}

	dx mkdir -p ${ODIR}

	BCF=/mnt/project/UKB_PHASING_EXOME_ARRAY/step0_merge/chr${CHR}/UKB.chr${CHR}.exome_array.WO_parents.sorted.bcf		
	OUT=UKB_chr${CHR}.exome_array.WO_parents.shapeit5.common_${CUT}.bcf
	LOG=UKB_chr${CHR}.exome_array.WO_parents.shapeit5.common_${CUT}.log
	TIM=UKB_chr${CHR}.exome_array.WO_parents.shapeit5.common_${CUT}.time
	
	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_220824_c85eb0c.tar.gz" --folder="./$ODIR/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase_common_static --input $BCF --map $MAP --output $OUT --thread $THREADS --log $LOG --filter-maf $CUT --region chr${CHR} && bcftools index -f $OUT --threads $THREADS" --tag shapeit5_common --tag chr$CHR --instance-type mem2_ssd1_v2_x64 --priority normal --name shapeit5_common_FULL_chr${CHR} -y

	

done
