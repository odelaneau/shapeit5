#!bin/bash

threads=8

for CHR in {1..22}; do


	mkdir -p TMP
	dx ls Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants,\ pVCF\ format\ -\ interim\ 450k\ release/ukb23148_c${CHR}_b*.gz > TMP/chr${CHR}.chunks


	ODIR=UKB_PHASING_EXOME_ARRAY/step0_merge/chr${CHR}/support
	for file in $(cat TMP/chr${CHR}.chunks); do
	
		chunk=$(echo $file | cut -d '_' -f 3 | cut -d 'b' -f 2)
		exome_overlapping_samples_chunk=ukb23148_c${CHR}_v1.b${chunk}.overlap_array.bcf
		overlapping_samples="/mnt/project/UKB_PHASING_EXOME_ARRAY/step0_merge/chr${CHR}/support/chr${CHR}.overlapping_samples"
		input="/mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants,\ pVCF\ format\ -\ interim\ 450k\ release/${file}"
		BCF="UKB.chr${CHR}.b${chunk}.exome_array.split_multiallelic.bcf"


		dx run app-swiss-army-knife -icmd="bcftools view --threads ${threads} -S ${overlapping_samples} -Ob -o ${exome_overlapping_samples_chunk} ${input} && bcftools index ${exome_overlapping_samples_chunk} && bcftools annotate -x ^FORMAT/GT ${exome_overlapping_samples_chunk} | bcftools annotate -x FILTER/MONOALLELIC | bcftools annotate -x ^INFO/AC,INFO/AN | bcftools reheader --threads ${threads} -h /mnt/project/UKB_PHASING_EXOME_ARRAY/step0_merge/chr${CHR}/support/array_header.chr${CHR}.vcf.gz > ukb23148_c${CHR}_v1.b${chunk}.overlap_array.TAGS.vcf && bgzip -f ukb23148_c${CHR}_v1.b${chunk}.overlap_array.TAGS.vcf && tabix -p vcf -f ukb23148_c${CHR}_v1.b${chunk}.overlap_array.TAGS.vcf.gz && bcftools norm -m -any --threads ${threads} -Ob -o tmp.bcf ukb23148_c${CHR}_v1.b${chunk}.overlap_array.TAGS.vcf.gz && bcftools index tmp.bcf && bcftools annotate --threads ${threads} --set-id '%CHROM\_%POS\_%REF\_%ALT' -Ob -o tmp2.bcf tmp.bcf && bcftools index tmp2.bcf && bcftools view -e 'ID=@/mnt/project/UKB_PHASING_EXOME_ARRAY/step0_merge/chr${CHR}/support/chr${CHR}.array_snps_kept.txt' -Ob -o ${BCF} tmp2.bcf && bcftools index ${BCF} && rm tmp* && rm ukb23148_c${CHR}_v1.b${chunk}.overlap_array.*" --tag chr${CHR} --tag chunk_${chunk} --instance-type mem1_ssd1_v2_x8 --folder="./${ODIR}" --name part1_process_exome_chr${CHR}_b${chunk} --priority normal -y









	done
done


