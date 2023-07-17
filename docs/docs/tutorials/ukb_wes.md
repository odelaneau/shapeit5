---
layout: default
title: UK Biobank WES data
nav_order: 4
parent: Tutorials
---
# UK Biobank WES data
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

**For any question on this pipeline, please contact [Robin J. Hofmeister](https://github.com/RJHFMSTR) and [Olivier Delaneau](https://github.com/odelaneau).** 
 
---

## Rationale
Similar to the WGS data, the WES data can be phased in two steps, separating common and rare markers. However, the specificity of the WES phasing is that we first merged it with the SNP array data to increase the density of common markers, in particular in between genes.
<br>

## Phasing the WES data
<br>
### Set up your environment
To structure the outputs of the analysis and to simplify the understanding of each step of this tutorial, we first create output folders as follows. You can choose the change the name of these folders but you will have to change the code accordingly.
<div class="code-example" markdown="1">
```bash
dx mkdir -p Phasing/PhasingWES/step0_merge/support/
dx mkdir -p Phasing/PhasingWES/step1_phase_common/
dx mkdir -p Phasing/PhasingWES/step2_phase_rare/chunks/
```
</div>
<br>

### Merging WES and SNP array datas

To merge WES and SNP array data, we proceed in several steps as follows.

Important note: This uses some files produced as part of the UK Biobank SNP array tutorial, such as the SNP QC list (that should be located here after following the tutorial Phasing/PhasingSNParray/step1_dataqc/SNPlist.filtered.QC.txt)

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
	dx run app-swiss-army-knife -icmd="bcftools query -l ${exome_chunks1} > chr${CHR}.exome_samples && bcftools query -l ${array} | cut -d '_' -f 1 > chr${CHR}.array_samples && cat chr${CHR}.exome_samples chr${CHR}.array_samples > chr${CHR}.samples && sort chr${CHR}.samples | uniq -c | gawk '\$1==2{print \$2}' > chr${CHR}.overlapping_samples && bcftools reheader --threads ${threads} -s chr${CHR}.array_samples ${array} | bcftools view --threads ${threads} -S chr${CHR}.overlapping_samples -Oz -o ${array_overlapping_samples} && bcftools index -c --threads ${threads} ${array_overlapping_samples} && rm chr${CHR}.exome_samples && rm chr${CHR}.array_samples && rm chr${CHR}.samples && bcftools annotate -x ^INFO/AC,INFO/AN ${array_overlapping_samples} | bcftools view -i 'ID=@/mnt/project//Phasing/PhasingSNParray/step1_dataqc/SNPlist.filtered.QC.txt' | bcftools annotate --threads ${threads} --set-id '%CHROM\_%POS\_%REF\_%ALT' -Oz -o full_c${CHR}_b0_v2.b38.sorted.overlap_exome.TAGS.vcf.gz && tabix -p vcf -f full_c${CHR}_b0_v2.b38.sorted.overlap_exome.TAGS.vcf.gz && rm full_c${CHR}_b0_v2.b38.sorted.overlap_exome.vcf.gz* && bcftools view -r chr20:1-10 -Oz -o array_header.chr${CHR}.vcf.gz full_c${CHR}_b0_v2.b38.sorted.overlap_exome.TAGS.vcf.gz && bcftools index array_header.chr${CHR}.vcf.gz && bcftools query -f '%ID\n' full_c${CHR}_b0_v2.b38.sorted.overlap_exome.TAGS.vcf.gz > chr${CHR}.array_snps_kept.txt" --tag chr${CHR} --tag overlap_samples --instance-type mem1_ssd1_v2_x16 --folder="./Phasing/PhasingWES/step0_merge/support/" --name overlap_indiv_chr${CHR}_SNParray --priority normal -y
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

<br>





### Phasing
SHAPEIT5 phases common variants using the SHAPEIT5_phase_common tool. As an input, **SHAPEIT5_phase_common** requires an unphased file (with AC and AN tags), and automatically sub-sets the file to the desired MAF (e.g. 0.001). There are different strategies to phase common variants. Conversly to the WGS phasing, the phasing of common variants for the WES data can be performed across entire chromosomes. However, the phasing of rare variants will be performed in chunks.
<br>

#### Phasing common
We phase common variants using **SHAPEIT5_phase_common** across entire chromosomes. This phasing is then used as a scaffold to phase rare variants in chunks.


**IMPORTANT**: in the following code make sure to change the shapeit5 docker image name (here `shapeit5_beta.tar.gz`) to the latest version that you've downloaded [here](https://odelaneau.github.io/shapeit5/docs/installation/docker)


<div class="code-example" markdown="1">
```bash
# step1. Download map files
dx mkdir -p data/shapeit_maps/
wget https://github.com/odelaneau/shapeit5/raw/main/maps/genetic_maps.b38.tar.gz
tar -xvzf genetic_maps.b38.tar.gz
dx upload *.b38.gmap.gz --path="data/shapeit_maps/"

# step2. Phasing
for CHR in {1..22}; do
	MAP=/mnt/project/data/shapeit_maps/chr${CHR}.b38.gmap.gz
	CUT=0.001
	THREADS=64
	BCF=/mnt/project/Phasing/PhasingWES/step0_merge/UKB.chr${CHR}.exome_array.full.sorted.bcf		
	OUT=UKB_chr${CHR}.exome_array.shapeit5.common_${CUT}.bcf
	LOG=UKB_chr${CHR}.exome_array.shapeit5.common_${CUT}.log
	TIM=UKB_chr${CHR}.exome_array.shapeit5.common_${CUT}.time
	dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_beta.tar.gz" --folder="./Phasing/PhasingWES/step1_phase_common/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase_common_static --input $BCF --map $MAP --output $OUT --thread $THREADS --log $LOG --filter-maf $CUT --region chr${CHR} && bcftools index -f $OUT --threads $THREADS" --instance-type mem2_ssd1_v2_x64 --priority normal --name shapeit5_common_chr${CHR} -y
done
```
</div>

The full list of options for the **SHAPEIT5_phase_common** command can be found [here](https://odelaneau.github.io/shapeit5/docs/documentation/phase_common/).


<br>
#### Phasing rare

We phase rare variants using **SHAPEIT5_phase_rare** for small regions (i.e chunks), that we then concatenate to resolve the phasing for the entire chromosome. For this, we use as a scaffold the phasing of common variants (previous step)


<div class="code-example" markdown="1">
```bash
# step1. Download chunk file
(available soon)

# step2. Phasing 
for CHR in {1..22}; do
	MAP=/mnt/project/data/shapeit_maps/chr${CHR}.b38.gmap.gz
	CHUNKS=
	CUT=0.001
	THREADS=32
	while read LINE; do
			SCAFFOLD_REG=$(echo $LINE | awk '{ print $3; }')
			SCAFFOLD_REG_START=$(echo ${SCAFFOLD_REG} | cut -d":" -f 2 | cut -d"-" -f1)
			SCAFFOLD_REG_END=$(echo ${SCAFFOLD_REG} | cut -d":" -f 2 | cut -d"-" -f2)
			SCAFFOLD_REG_NAME=${CHR}_${SCAFFOLD_REG_START}_${SCAFFOLD_REG_END}
			INPUT_REG=$(echo $LINE | awk '{ print $4; }')
			INPUT_REG_START=$(echo ${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f1)
			INPUT_REG_END=$(echo ${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f2)
			INPUT_REG_NAME=${CHR}_${INPUT_REG_START}_${INPUT_REG_END}
			BCF=/mnt/project/Phasing/PhasingWES/step0_merge/UKB.chr${CHR}.exome_array.WO_parents.sorted.bcf
			SCAFFOLD=/mnt/project/Phasing/PhasingWES/step1_phase_common/UKB_chr${CHR}.exome_array.shapeit5.common_${CUT}.bcf
			OUT=UKB_chr${CHR}.exome_array.${INPUT_REG_NAME}.shapeit5.WO_parents.rares_${CUT}.bcf
			LOG=UKB_chr${CHR}.exome_array.${INPUT_REG_NAME}.shapeit5.WO_parents.rares_${CUT}.log
			TIM=UKB_chr${CHR}.exome_array.${INPUT_REG_NAME}.shapeit5.WO_parents.rares_${CUT}.time
			dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_beta.tar.gz" --folder="./Phasing/PhasingWES/step2_phase_rare/chunks/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase_rare_static --input $BCF --scaffold $SCAFFOLD --map $MAP --output $OUT --log $LOG --scaffold-region $SCAFFOLD_REG --input-region $INPUT_REG --thread $THREADS && bcftools index -f $OUT --threads $THREADS" --instance-type mem2_ssd1_v2_x32 --priority normal --name shapeit5_rares_$INPUT_REG_NAME -y		
	done < $CHUNKS
done
```
</div>


To accurately phase the edges of our chunks, we included a buffer region, overlapping the previous and next chunks. Before merging all phasing chunks, we need to remove this buffer to avoid aving duplicated markers.

<div class="code-example" markdown="1">
```bash
for CHR in {1..22}; do
	MAP=/mnt/project/data/shapeit_maps/chr${CHR}.b38.gmap.gz
	CHUNKS=
	CUT=0.001
	THREADS=16
	while read LINE; do
			INPUT_REG=$(echo $LINE | awk '{ print $4; }')
			INPUT_REG_START=$(echo ${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f1)
			INPUT_REG_END=$(echo ${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f2)
			INPUT_REG_NAME=${CHR}_${INPUT_REG_START}_${INPUT_REG_END}
			INPUT_NBR=$(echo $LINE | awk '{ print $1; }')
			IN=/mnt/project//Phasing/PhasingWES/step2_phase_rare/chunks/UKB_chr${CHR}.exome_array.${INPUT_REG_NAME}.shapeit5.rares_${CUT}.bcf
			OUT=UKB_chr${CHR}.exome_array.${INPUT_NBR}.shapeit5.rares_${CUT}.bcf
			dx run app-swiss-army-knife --folder="./Phasing/PhasingWES/step2_phase_rare/chunks/" -icmd="bcftools view -r ${INPUT_REG} --threads ${THREADS} -Ob -o ${OUT} ${IN} && bcftools index ${OUT} --threads ${THREADS}" --instance-type mem1_ssd1_v2_x16 --priority normal --name trim_$INPUT_REG_NAME -y
	done < $CHUNKS
done
```
</div>



We can finally merge all our phasing chunks using the command ""bcftools concat**.


<div class="code-example" markdown="1">
```bash
for CHR in {1..22}; do
	dx run app-swiss-army-knife --folder="./Phasing/PhasingWES/step2_phase_rare/" -icmd="bcftools concat -n /mnt/project/Phasing/PhasingWES/step2_phase_rare/chunks/UKB_chr${CHR}.exome_array.*.shapeit5.rares_*.bcf" -Ob -o UKB_chr${CHR}.phased.bcf && bcftools index UKB_chr${CHR}.phased.bcf" --instance-type mem1_ssd1_v2_x16 --priority normal --name concat_chr${CHR} -y
done
```
</div>




## Validation of your phasing
You can validate the quality of the haplotypes using the **SHAPEIT5_switch** tool. For this you will need parent-offspring duos or trios stored in a three-columns file (here called `family.ped`) following this format:

<div class="code-example" markdown="1">
- family.ped:      `offspring_id`    `parent1_id`    `parent2_id`

</div>


To validate the phasing using family data, the phasing must be performed by excluding parental genomes, so that offsprings are phased regardless of their parental genomes. This can be done using the **bcftools view** command, with as input the merged SNP array + WES data (located here `/mnt/project/Phasing/PhasingWES/step0_merge/UKB.chr${CHR}.exome_array.full.sorted.bcf`).

<div class="code-example" markdown="1">
```bash
dx mkdir -p Phasing/PhasingWES/benchmark/
for CHR in {1..22}; do
IN=/mnt/project/Phasing/PhasingWES/step0_merge/UKB.chr${CHR}.exome_array.full.sorted.bcf
OUT=benchmark_UKB.chr${CHR}.exome_array.full.sorted.bcf
dx run app-swiss-army-knife --folder "/Phasing/PhasingWES/benchmark/" -icmd="bcftools view --threads 16 -S ^parents.txt -Ob -o ${OUT} ${IN} && bcftools index ${OUT} --threads 16" --instance-type mem1_ssd1_v2_x16 --priority normal --name benchmark_wes_chr${CHR} -y
done
```
</div>

After excluding parental genomes with the above command, proceed with the normal phasing procedure as detailed above.

Let's consider that you performed the above steps of phasing using the input data exluding parental genomes, which produced a phased output file that you named `benchmark_UKB.chr${CHR}.exome_array.full.sorted.phased.bcf`). You can validate you phasing using the following command:


<div class="code-example" markdown="1">
```bash
for CHR in 20; do
dx run app-swiss-army-knife --folder "/Phasing/PhasingWES/benchmark/" -iimage_file="/docker/shapeit5_beta.tar.gz" -icmd="SHAPEIT5_switch --validation /mnt/project/Phasing/PhasingWES/step0_merge/UKB.chr${CHR}.exome_array.full.sorted.bcf --estimation benchmark_UKB.chr${CHR}.exome_array.full.sorted.phased.bcf --region ${CHR} --output benchmark_wes_chr${CHR}"  --instance-type mem2_ssd1_v2_x16 --priority normal --name benchmark_wes_chr${CHR} -y
done
```
</div>



The full list of options for the **SHAPEIT5_switch** command can be found [here](https://odelaneau.github.io/shapeit5/docs/documentation/switch/).















