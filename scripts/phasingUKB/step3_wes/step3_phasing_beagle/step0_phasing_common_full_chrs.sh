#!/bin/bash




threads=96


for CHR in 2; do

	ODIR=UKB_PHASING_EXOME_ARRAY/step0_merge/chr${CHR}/
	BCF=/mnt/project/UKB_PHASING_EXOME_ARRAY/step0_merge/chr${CHR}/UKB.chr${CHR}.exome_array.WO_parents.sorted.bcf	
	OUT=UKB.chr${CHR}.exome_array.WO_parents.sorted.vcf.gz
	
	#dx run app-swiss-army-knife --folder="./$ODIR/" -icmd="bcftools view --threads $threads} -Oz -o ${OUT} ${BCF} && bcftools index -f ${OUT} --threads ${threads}" --tag chr$CHR --instance-type mem1_ssd1_v2_x36 --priority normal --name beage_convert_chr${CHR} -y
	
	
	VCF=/mnt/project/UKB_PHASING_EXOME_ARRAY/step0_merge/chr${CHR}/${OUT}
	BGL=/mnt/project/docker/beagle.19Apr22.7c0.jar
	MAP=/mnt/project/data/plink_maps/plink.prefix.chr${CHR}.GRCh38.map	
	CHUNKS=../step1_chunking/chunks.chr${CHR}.txt
	CUT=0.01
	THREADS=64
	#ODIR=UKB_PHASING_EXOME_ARRAY/step5_phasing_beagle_common/chr${CHR}

	dx mkdir -p ${ODIR}

	
	
	OUT=UKB_chr${CHR}.exome_array.WO_parents.beagle.common.vcf.gz
	LOG=UKB_chr${CHR}.exome_array.WO_parents.beagle.common.log
	TIM=UKB_chr${CHR}.exome_array.WO_parents.beagle.common.time
	
	
	dx run app-swiss-army-knife --folder="./$ODIR/" -icmd="/usr/bin/time -vo $TIM java -Xmx256G -jar $BGL gt=$VCF map=$MAP out=$OUT nthreads=${threads} chrom=chr${CHR} && bcftools index -f $OUT\.vcf.gz --threads 32" --tag benchWGS --tag beagle5.4 --tag chr20 --instance-type mem2_ssd1_v2_x96 --priority normal --name beagle_chr${CHR} -y
	
	

done
