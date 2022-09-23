#!bin/bash

threads=36


for CHR in {1..22}; do
	ODIR=UKB_PHASING_EXOME_ARRAY/step0_merge/chr${CHR}




	array_vcf="/mnt/project/UKB_PHASING_EXOME_ARRAY/step0_merge/chr${CHR}/support/full_c${CHR}_b0_v2.b38.sorted.overlap_exome.TAGS.vcf.gz"
	array_bcf="full_c${CHR}_b0_v2.b38.sorted.overlap_exome.TAGS.bcf"
	chunks="/mnt/project/UKB_PHASING_EXOME_ARRAY/step0_merge/chr${CHR}/support/UKB.chr${CHR}.b*.exome_array.split_multiallelic.bcf"
	
	merged=UKB.chr${CHR}.exome_array.full.sorted.bcf


	tmp1=tmp1.bcf
	tmp2=tmp2.bcf



	#dx run app-swiss-army-knife -icmd="bcftools view --threads ${threads} -Ob -o ${array_bcf} ${array_vcf} && bcftools index ${array_bcf} --threads ${threads} && bcftools concat --naive-force --threads ${threads} ${array_bcf} ${chunks} -Ob -o ${tmp1} && bcftools sort -Ob -o ${tmp2} ${tmp1} && bcftools index ${tmp2} --threads ${threads} && rm ${tmp1} && bcftools view --threads ${threads} -i 'F_MISSING < 0.1' -Ob -o ${merged} ${tmp2} && bcftools index ${merged} --threads ${threads} && rm ${tmp2}*" --tag chr${CHR} --tag WESarray_merge --instance-type mem1_hdd1_v2_x36 --folder="./${ODIR}" --name WESarray_merge_chr${CHR} --priority low -y



	# remove parents
	parents=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/samples.parents.txt
	IN=/mnt/project/UKB_PHASING_EXOME_ARRAY/step0_merge/chr${CHR}/UKB.chr${CHR}.exome_array.full.sorted.bcf
	OUT=UKB.chr${CHR}.exome_array.WO_parents.sorted.bcf

	dx run app-swiss-army-knife -icmd="bcftools view -S^${parents} --force-samples --threads ${threads} -Ob -o ${OUT} ${IN} && bcftools index ${OUT} --threads ${threads}" --tag chr${CHR} --instance-type mem1_hdd1_v2_x36 --folder="./${ODIR}" --name WESremoveparents_chr${CHR} --priority low -y


done


