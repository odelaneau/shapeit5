---
layout: default
title: UK Biobank WES data
nav_order: 4
parent: Tutorials
---
# UK Biobank WES data
{: .no_toc }

{: .warning }
Website under construction: content not available yet!

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---


## Rationale
Similar to the WGS data, the WES data can be phased in two steps, separating common and rare markers. However, the specificity of the WES phasing is that we first merged it with the SNP array data to increase the density of common markers, in particular in between genes.

## Phasing the WES data

dx mkdir -p Phasing/PhasingWES/step0_merge/support/

### Merging WES and SNP array datas

To merge WES and SNP array data, we proceed in several steps as follows:

**1. SNP array lifting over.**

To merge the SNP array with the WES, we first quality control the data and lift it over to hg38. This is described in the [*UK Biobank SNP array data* tutorial](https://odelaneau.github.io/shapeit5/docs/tutorials/ukb_snp_array/).

**2. Subsetting overlapping individuals.**

In that step we keep only individuals listed in both the SNP array and the WES data

<div class="code-example" markdown="1">
```bash
# step1. Get overlapping individuals and subset the SNP array data
for CHR in {1..22}; do
	exome_chunks1="/mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants,\ pVCF\ format\ -\ interim\ 450k\ release/ukb23148_c${CHR}_b1_v1.vcf.gz"
	array="/Phasing/PhasingSNParray/step4_liftover/full_c$CHR\_b0_v2.b38.sorted.vcf.gz"
	array_overlapping_samples=full_c${CHR}_b0_v2.b38.sorted.overlap_exome.vcf.gz
	threads=16
	dx run app-swiss-army-knife -icmd="bcftools query -l ${exome_chunks1} > chr${CHR}.exome_samples && bcftools query -l ${array} | cut -d '_' -f 1 > chr${CHR}.array_samples && cat chr${CHR}.exome_samples chr${CHR}.array_samples > chr${CHR}.samples && sort chr${CHR}.samples | uniq -c | gawk '\$1==2{print \$2}' > chr${CHR}.overlapping_samples && bcftools reheader --threads ${threads} -s chr${CHR}.array_samples ${array} | bcftools view --threads ${threads} -S chr${CHR}.overlapping_samples -Oz -o ${array_overlapping_samples} && bcftools index -c --threads ${threads} ${array_overlapping_samples} && rm chr${CHR}.exome_samples && rm chr${CHR}.array_samples && rm chr${CHR}.samples && bcftools annotate -x ^INFO/AC,INFO/AN ${array_overlapping_samples} | bcftools view -i 'ID=@/mnt/project/UKB_PHASING_EXOME_ARRAY/SNPs_in_phasing.QC.txt' | bcftools annotate --threads ${threads} --set-id '%CHROM\_%POS\_%REF\_%ALT' -Oz -o full_c${CHR}_b0_v2.b38.sorted.overlap_exome.TAGS.vcf.gz && tabix -p vcf -f full_c${CHR}_b0_v2.b38.sorted.overlap_exome.TAGS.vcf.gz && rm full_c${CHR}_b0_v2.b38.sorted.overlap_exome.vcf.gz* && bcftools view -r chr20:1-10 -Oz -o array_header.chr${CHR}.vcf.gz full_c${CHR}_b0_v2.b38.sorted.overlap_exome.TAGS.vcf.gz && bcftools index array_header.chr${CHR}.vcf.gz && bcftools query -f '%ID\n' full_c${CHR}_b0_v2.b38.sorted.overlap_exome.TAGS.vcf.gz > chr${CHR}.array_snps_kept.txt" --tag chr${CHR} --tag overlap_samples --instance-type mem1_ssd1_v2_x16 --folder="./Phasing/PhasingWES/step0_merge/support/" --name overlap_indiv_chr${CHR}_SNParray --priority normal -y
done


# step2. Subset the WES data, split multi-allelic sites, remove markers present in both the WES and the SNP array data
for CHR in {1..22}; do
	threads=8
	mkdir -p TMP
	dx ls Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants,\ pVCF\ format\ -\ interim\ 450k\ release/ukb23148_c${CHR}_b*.gz > TMP/chr${CHR}.chunks

	for file in $(cat TMP/chr${CHR}.chunks); do
		chunk=$(echo $file | cut -d '_' -f 3 | cut -d 'b' -f 2)
		exome_overlapping_samples_chunk=ukb23148_c${CHR}_v1.b${chunk}.overlap_array.bcf
		overlapping_samples="/mnt/project/Phasing/PhasingWES/step0_merge/support/chr${CHR}.overlapping_samples"
		input="/mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants,\ pVCF\ format\ -\ interim\ 450k\ release/${file}"
		BCF="UKB.chr${CHR}.b${chunk}.exome_array.split_multiallelic.bcf"
		dx run app-swiss-army-knife -icmd="bcftools view --threads ${threads} -S ${overlapping_samples} -Ob -o ${exome_overlapping_samples_chunk} ${input} && bcftools index ${exome_overlapping_samples_chunk} && bcftools annotate -x ^FORMAT/GT ${exome_overlapping_samples_chunk} | bcftools annotate -x FILTER/MONOALLELIC | bcftools annotate -x ^INFO/AC,INFO/AN | bcftools reheader --threads ${threads} -h /mnt/project/UKB_PHASING_EXOME_ARRAY/step0_merge/chr${CHR}/support/array_header.chr${CHR}.vcf.gz > ukb23148_c${CHR}_v1.b${chunk}.overlap_array.TAGS.vcf && bgzip -f ukb23148_c${CHR}_v1.b${chunk}.overlap_array.TAGS.vcf && tabix -p vcf -f ukb23148_c${CHR}_v1.b${chunk}.overlap_array.TAGS.vcf.gz && bcftools norm -m -any --threads ${threads} -Ob -o tmp.bcf ukb23148_c${CHR}_v1.b${chunk}.overlap_array.TAGS.vcf.gz && bcftools index tmp.bcf && bcftools annotate --threads ${threads} --set-id '%CHROM\_%POS\_%REF\_%ALT' -Ob -o tmp2.bcf tmp.bcf && bcftools index tmp2.bcf && bcftools view -e 'ID=@/mnt/project/UKB_PHASING_EXOME_ARRAY/step0_merge/chr${CHR}/support/chr${CHR}.array_snps_kept.txt' -Ob -o ${BCF} tmp2.bcf && bcftools index ${BCF} && rm tmp* && rm ukb23148_c${CHR}_v1.b${chunk}.overlap_array.*" --tag chr${CHR} --tag chunk_${chunk} --instance-type mem1_ssd1_v2_x8 --folder="./Phasing/PhasingWES/step0_merge/support/" --name overlap_indiv_chr${CHR}_b${chunk}_WES --priority normal -y
	done
done
```
</div>


**3. Merge.**


<div class="code-example" markdown="1">
```bash
for CHR in {1..22}; do
	threads=36
	array_vcf="/mnt/project/Phasing/PhasingWES/step0_merge/support/full_c${CHR}_b0_v2.b38.sorted.overlap_exome.TAGS.vcf.gz"
	array_bcf="full_c${CHR}_b0_v2.b38.sorted.overlap_exome.TAGS.bcf"
	chunks="/mnt/project/Phasing/PhasingWES/step0_merge/support/UKB.chr${CHR}.b*.exome_array.split_multiallelic.bcf"
	merged=UKB.chr${CHR}.exome_array.full.sorted.bcf
	dx run app-swiss-army-knife -icmd="bcftools view --threads ${threads} -Ob -o ${array_bcf} ${array_vcf} && bcftools index ${array_bcf} --threads ${threads} && bcftools concat --naive-force --threads ${threads} ${array_bcf} ${chunks} -Ob -o tmp1.bcf && bcftools sort -Ob -o tmp2.bcf tmp1.bcf && bcftools index tmp2.bcf --threads ${threads} && rm tmp1.bcf* && bcftools view --threads ${threads} -i 'F_MISSING < 0.1' -Ob -o ${merged} tmp2.bcf && bcftools index ${merged} --threads ${threads} && rm tmp2.bcf*" --instance-type mem1_hdd1_v2_x36 --folder="./Phasing/PhasingWES/step0_merge" --name WESarray_merge_chr${CHR} --priority low -y

```
</div>























