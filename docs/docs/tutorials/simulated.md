---
layout: default
title: Simulated WGS data
nav_order: 1
parent: Tutorials
---
# Simulated WGS data
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Rationale
SHAPEIT5 is a two-step approach that treats each chromosome independently and works as follows:

1. Phase common variants (MAF >= 0.1%) of a chromosome using **SHAPEIT5_phase_common**. This can be done as a single job for SNP array or WES data, but it might be necessary to split the chromosome into large chunks (e.g. 20 cM) in the case of WGS data.

2. Ligate the phased common variants (MAF >= 0.1%) of a chromosome using **SHAPEIT5_ligate**, only if chunking was performed in step 1. The ligation step is computationally light and uses variants in the intersection of the chunks to provide chromosome-wide haplotypes. The result of this step (or the previous step if no chunking was used), is used as a haplotype scaffold for the next step.

3. Phase rare variants (MAF < 0.1%) of a chromosome using **SHAPEIT5_rare**. To do this, we use the haplotype scaffold generated in step 2 (or 1) and we proceed in relatively small chunks (e.g. 5Mb) to run relatively fast job in parallel. At the end of this step, we have several fully phased chunks across the chromosome.

4. Concatenate the phased chunks generated in step 3 using **bcftools concat -n**. As in the previous step haplotypes have been phased onto a haplotype scaffold, there is no need to ligate the chunks, and the files can be concatenated without decompression and recompression. This makes this step almost instantaneous, even for large cohorts.

## Phasing of simulated data
For this tutorial, we work on a dense 10Mb region generated with msprime for 100,000 European samples. The data is perfectly phased, therefore we will be able to perform haplotype phasing on a large cohort and also to verify the accuracy of the phasing process. **IMPORTANT**: as we are using simulated data here, we do not provide a recombination map in the following steps. However, this is **very important** when using the method on real data.

### Phasing common variants
SHAPEIT5 phases common variants using the **SHAPEIT5_phase_common** tool. As an input, SHAPEIT5_phase_common requires an unphased file (with AC and AN tags), and automatically sub-sets the file to the desired MAF (e.g. 0.001).

There are different strategies to phase common variants. The first, is to phase the whole chromosome in a single job. This is feasible for SNP array data and WES data, but it is not optimal for WGS data. Therefore we recommend to chunk the chromosome into large chunks (e.g. 20 cM) if using large WGS data. In the following we see how to impute the simulated region using both strategies: the chunking performed in this tutorial is under optimal for real data, such as the UK Biobank dataset. For that please refer to the tutorial using real data.


Before starting this tutorial, be sure to clone the SHAPEIT5 github and compile SHAPEIT5 (documentation [here](https://odelaneau.github.io/shapeit5/docs/installation/build_from_source/compile_shapeit5), and to navigate in the main shapeit5 folder. In that folder, you should see at the following folders:

<div class="code-example" markdown="1">
```bash
$ ls -l

common/
docker/
ligate/
maps/
phase_common/
phase_rare/
switch/
test/
```
</div>


The simulated data that we use in this tutorial are located in `test/10k/`.

#### Phasing all common variants in one job
SHAPEIT5 can phase common variants (MAF >= 0.001) for the the whole 10Mb region using the following command:

<div class="code-example" markdown="1">
```bash
../phase_common/bin/SHAPEIT5_phase_common --input 10k/msprime.nodup.bcf --filter-maf 0.001  --output 10k/msprime.common.phased.bcf --region 1 --thread 8
bcftools index 10k/msprime.common.phased.bcf --threads 8
```
</div>

As this is the entire region of the simulation, **no ligation step is required**, and therefore in this scenario we can directly validate our haplotypes at common variants and proceed with rare variant phasing.

#### Phasing common variants in chunks
More realistically, when using WGS data on large sample size, it is good practice to run **SHAPEIT5_phase_common** in different large regions of the chromosomes (e.g. 20cM). In the following, we perform phasing in two overlapping 6Mb windows. The intersection of the two chunks is therefore 2Mb, that is large enough to have a good amount of heterozygous sites for the ligation step.

<div class="code-example" markdown="1">
```bash
#first chunk: 1:1-6000000
../phase_common/bin/SHAPEIT5_phase_common --input 10k/msprime.nodup.bcf --filter-maf 0.001  --output 10k/msprime.common.phased_chr1_1_6000000.bcf --region 1:1-6000000 --thread 8
bcftools index 10k/msprime.common.phased_chr1_1_6000000.bcf --threads 8

#second chunk: 1:4000001-10000000
../phase_common/bin/SHAPEIT5_phase_common --input 10k/msprime.nodup.bcf --filter-maf 0.001  --output 10k/msprime.common.phased_chr1_4000001_10000000.bcf --region 1:4000001-10000000 --thread 8
bcftools index 10k/msprime.common.phased_chr1_4000001_10000000.bcf --threads 8
```
</div>

One advantage of this approach is that these two chunks can run in parallel and they require less computational resources than the single job. However, as we wish to perform phasing at rare variants using a haplotype scaffold, we need to ligate these chunks to create a single file for the entire region.

#### Ligate haplotypes at common variants
Ligation  of the phased common variants is performed using the **SHAPEIT5_ligate** tool. The program requires an ordered list of the phased common variants. We recommend using appropriate naming of the files in the previous step, so that a command such as `ls -1v` can directly produce the list of files in the right order.

<div class="code-example" markdown="1">
```bash
#generate list for SHAPEIT5_ligate
ls -1v 10k/msprime.common.phased_chr1_*.bcf > list_ligate.txt
#ligate the chunks
../phase_ligate/bin/SHAPEIT5_ligate --input list_ligate.txt --output 10k/msprime.ligated.phased.bcf --thread 8 --index
```
</div>

At the end of this step, we have a region-wide haplotype scaffold.

### Validation of haplotypes at common variants
We can validate the quality of the haplotypes at common variants using the **SHAPEIT5_switch** tool. For this we use the known phased haplotypes as validation dataset, and the (ligated) region-wide file phased in the previous steps.

<div class="code-example" markdown="1">
```bash
../switch/bin/SHAPEIT5_switch --validation 10k/msprime.nodup.bcf --estimation 10k/msprime.common.phased.bcf --region 1 --output 10k/msprime.common.phased
```
</div>


### Phase rare variants small region
Having obtained the haplotype scaffold at common variants in the previous step, we can now phase rare variants using the **SHAPEIT5_rare** tool. We always do this in chunks, in order to maximise the parallelization of our jobs (e.g. in 5Mb chunks). For this simulated data, we use ~1Mb chunks, with a buffer of 0.5Mb at each side. The program requires two input files, one is the haplotype scaffold generated in the previous steps, and other input require is a file containing the whole unphased region, the same file provided for **SHAPEIT5_phase_common**: SHAPEIT5 automatically discards common variants for this file, duplicated in the phased scaffold.

<div class="code-example" markdown="1">
```bash
#chunk1
./phase_rare/bin/SHAPEIT5_phase_rare --input-plain 10k/msprime.nodup.bcf --scaffold 10k/msprime.common.phased.bcf --output 10k/msprime.rare.chr1_1_1500000.bcf --scaffold-region 1:1-2000000 --input-region 1:1-1500000 --thread 8
bcftools index 10k/msprime.rare.chr1_1_1500000.bcf --threads 8
#chunk2
./phase_rare/bin/SHAPEIT5_phase_rare --input-plain 10k/msprime.nodup.bcf --scaffold 10k/msprime.common.phased.bcf --output 10k/msprime.rare.chr1_1500001_2500000.bcf --scaffold-region 1:1000000-3000000 --input-region 1:1500001-2500000 --thread 8
bcftools index 10k/msprime.rare.chr1_1500001_2500000.bcf --threads 8
#chunk3
./phase_rare/bin/SHAPEIT5_phase_rare --input-plain 10k/msprime.nodup.bcf --scaffold 10k/msprime.common.phased.bcf --output 10k/msprime.rare.chr1_2500001_3500000.bcf --scaffold-region 1:2000000-4000000 --input-region 1:2500001-3500000 --thread 8
bcftools index 10k/msprime.rare.chr1_2500001_3500000.bcf --threads 8
#chunk4
./phase_rare/bin/SHAPEIT5_phase_rare --input-plain 10k/msprime.nodup.bcf --scaffold 10k/msprime.common.phased.bcf --output 10k/msprime.rare.chr1_3500001_4500000.bcf --scaffold-region 1:3000000-5000000 --input-region 1:3500001-4500000 --thread 8
bcftools index 10k/msprime.rare.chr1_3500001_4500000.bcf --threads 8
#chunk5
./phase_rare/bin/SHAPEIT5_phase_rare --input-plain 10k/msprime.nodup.bcf --scaffold 10k/msprime.common.phased.bcf --output 10k/msprime.rare.chr1_4500001_5500000.bcf --scaffold-region 1:4000000-6000000 --input-region 1:4500001-5500000 --thread 8
bcftools index 10k/msprime.rare.chr1_4500001_5500000.bcf --threads 8
#chunk6
./phase_rare/bin/SHAPEIT5_phase_rare --input-plain 10k/msprime.nodup.bcf --scaffold 10k/msprime.common.phased.bcf --output 10k/msprime.rare.chr1_5500001_6500000.bcf --scaffold-region 1:5000000-7000000 --input-region 1:5500001-6500000 --thread 8
bcftools index 10k/msprime.rare.chr1_5500001_6500000.bcf --threads 8
#chunk7
./phase_rare/bin/SHAPEIT5_phase_rare --input-plain 10k/msprime.nodup.bcf --scaffold 10k/msprime.common.phased.bcf --output 10k/msprime.rare.chr1_6500001_7500000.bcf --scaffold-region 1:6000000-8000000 --input-region 1:6500001-7500000 --thread 8
bcftools index 10k/msprime.rare.chr1_6500001_7500000.bcf --threads 8
#chunk8
./phase_rare/bin/SHAPEIT5_phase_rare --input-plain 10k/msprime.nodup.bcf --scaffold 10k/msprime.common.phased.bcf --output 10k/msprime.rare.chr1_7500001_8500000.bcf --scaffold-region 1:7000000-9000000 --input-region 1:7500001-8500000 --thread 8
bcftools index 10k/msprime.rare.chr1_7500001_8500000.bcf --threads 8
#chunk9
./phase_rare/bin/SHAPEIT5_phase_rare --input-plain 10k/msprime.nodup.bcf --scaffold 10k/msprime.common.phased.bcf --output 10k/msprime.rare.chr1_8500001_1000000.bcf --scaffold-region 1:8000000-10000000 --input-region 1:8500001-1000000 --thread 8
bcftools index 10k/msprime.rare.chr1_8500001_1000000.bcf --threads 8
```
</div>

All these jobs can run in parallel. By default, **SHAPEIT5_phase_rare** discards buffer regions and therefore the output of a chunk does not overlap with others. Therefore, at the end of this step, we can simply concatenate the chunks in a straightforward way.

### Obtaining chromosome-wide phased data
As in the previous step haplotypes have been phased onto a haplotype scaffold, there is no need to ligate the chunks, and the files can be concatenated without decompression and recompression. This is performed very quickly with [bcftools concat \-\-naive](https://samtools.github.io/bcftools/bcftools.html#concat). We again recommend naming the files in an appropriate way, so that a command such as `ls -1v` can directly produce the list of files in the right order.

<div class="code-example" markdown="1">
```bash
#generate list for bcftools concat
ls -1v 10k/msprime.rare.chr1_*.bcf > list_concat.txt
#ligate the chunks
bcftools concat --naive -f list_concat.txt -o msprime.rare.phased.chr1.bcf --threads 8
bcftools index -f msprime.rare.phased.chr1.bcf --threads 8
```
</div>

### Validate phasing at rare variants
AS a final step, we can verify the quality of our haplotypes at rare variants using the **SHAPEIT5_switch** tool. Similarly to before we used the known haplotypes as validation and the phased rare variants. For WGS data, in order to speed-up the SWE calculations, you can run different **SHAPEIT5_switch** jobs in parallel, one per chunk. This is an approximation of the SWE at the entire region level, but in our tests the results are pretty consistent.

<div class="code-example" markdown="1">
```bash
../switch/bin/SHAPEIT5_switch --validation 10k/msprime.nodup.bcf --estimation 10k/msprime.rare.chr1_1_1500000.bcf --region 1 --output 10k/msprime.rare.chr1_1_1500000
```
</div>

