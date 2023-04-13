#!bin/bash

threads=16

for CHR in {1..22}; do

	ODIR=UKB_PHASING_EXOME_ARRAY/step0_merge/chr${CHR}/support

	exome_chunks1="/mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants,\ pVCF\ format\ -\ interim\ 450k\ release/ukb23148_c${CHR}_b1_v1.vcf.gz"
	array="/mnt/project/Phasing/PhasingSNParray/step4_liftover/full_c${CHR}_b0_v2.b38.sorted.vcf.gz"
	array_overlapping_samples=full_c${CHR}_b0_v2.b38.sorted.overlap_exome.vcf.gz
	
	dx run app-swiss-army-knife -icmd="bcftools query -l ${exome_chunks1} > chr${CHR}.exome_samples && bcftools query -l ${array} | cut -d '_' -f 1 > chr${CHR}.array_samples && cat chr${CHR}.exome_samples chr${CHR}.array_samples > chr${CHR}.samples && sort chr${CHR}.samples | uniq -c | gawk '\$1==2{print \$2}' > chr${CHR}.overlapping_samples && bcftools reheader --threads ${threads} -s chr${CHR}.array_samples ${array} | bcftools view --threads ${threads} -S chr${CHR}.overlapping_samples -Oz -o ${array_overlapping_samples} && bcftools index -c --threads ${threads} ${array_overlapping_samples} && rm chr${CHR}.exome_samples && rm chr${CHR}.array_samples && rm chr${CHR}.samples && bcftools annotate -x ^INFO/AC,INFO/AN ${array_overlapping_samples} | bcftools view -i 'ID=@/mnt/project/UKB_PHASING_EXOME_ARRAY/SNPs_in_phasing.QC.txt' | bcftools annotate --threads ${threads} --set-id '%CHROM\_%POS\_%REF\_%ALT' -Oz -o full_c${CHR}_b0_v2.b38.sorted.overlap_exome.TAGS.vcf.gz && tabix -p vcf -f full_c${CHR}_b0_v2.b38.sorted.overlap_exome.TAGS.vcf.gz && rm full_c${CHR}_b0_v2.b38.sorted.overlap_exome.vcf.gz* && bcftools view -r chr20:1-10 -Oz -o array_header.chr${CHR}.vcf.gz full_c${CHR}_b0_v2.b38.sorted.overlap_exome.TAGS.vcf.gz && bcftools index array_header.chr${CHR}.vcf.gz && bcftools query -f '%ID\n' full_c${CHR}_b0_v2.b38.sorted.overlap_exome.TAGS.vcf.gz > chr${CHR}.array_snps_kept.txt" --tag chr${CHR} --tag overlap_samples --instance-type mem1_ssd1_v2_x16 --folder="./${ODIR}" --name part0_process_array_chr${CHR} --priority normal -y

done


