#!/bin/bash



for CHR in 20 21 22; do


	CUT=0.001
	ODIR=UKB_PHASING_EXOME_ARRAY/step4_validation/Rares/chr${CHR}

	dx mkdir -p ${ODIR}

	
	

	INPUT=/mnt/project/UKB_PHASING_EXOME_ARRAY/step3_phasing_rares/UKB_chr${CHR}.exome_array.full.shapeit5.WO_parents.rares_${CUT}.bcf
	families=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/samples.families.txt
	BCF_estimation=UKB_chr${CHR}.trios_kids.bcf
	threads=36	

	#dx run app-swiss-army-knife -icmd="cat ${families} | grep -v 'NA' > trios.families.txt && cat trios.families.txt | cut -f 1 > trios.kids.txt && bcftools view --threads ${threads} --force-samples -S trios.kids.txt -Ob -o ${BCF_estimation} ${INPUT} && bcftools index ${BCF_estimation} --threads ${threads}" --tag chr${CHR} --tag bench_${chunk} --instance-type mem1_ssd1_v2_x36 --folder="./${ODIR}" --name bench1_${CHR} --priority normal -y



		



	#BCF=/mnt/project/UKB_PHASING_EXOME_ARRAY/step4_validation/Rares/chr${CHR}/${BCF_estimation}
	BCF=/mnt/project/UKB_PHASING_EXOME_ARRAY/step3_phasing_rares/UKB_chr${CHR}.exome_array.full.shapeit5.WO_parents.rares_0.001.bcf
	VAL=/mnt/project/UKB_PHASING_EXOME_ARRAY/OLD/step4_validation/Scaffold/chr${CHR}/validation_chr${CHR}.bcf
	FRQ=/mnt/project/UKB_PHASING_EXOME_ARRAY/step1_chunks/chr${CHR}/UKB.chr${CHR}.exome_array.sorted.sites.bcf
	PED=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/samples.families.txt
	OUT=SER_v2.chr${CHR}.cut${CUT}.fq
	LOG=SER_v2.chr${CHR}.cut${CUT}.fq.log
	threads=8

	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_2022-08-30_12bfbd7.tar.gz" --folder="./${ODIR}/" -icmd="SHAPEIT5_switch_static --validation $VAL --estimation $BCF --frequency $FRQ --pedigree $PED --region chr$CHR --output $OUT --log $LOG --thread $threads" --tag benchWGS --tag switchSHP4 --instance-type mem2_ssd2_v2_x8 --priority low --name benchWGS_chr${CHR} -y


		

done







