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

**For any question on imputation, please contact [Simone Rubinacci](https://srubinacci.github.io/).** 

---

## Rationale
SHAPEIT5 is a two-step approach that treats each chromosome independently and works as follows:

1. Phase common variants (MAF >= 0.1%) of a chromosome using **phase_common**. This can be done as a single job for SNP array, but it might be necessary to split the chromosome into large chunks (e.g. 20 cM) in the case of WGS data.

2. Ligate the phased common variants (MAF >= 0.1%) of a chromosome using **ligate**, only if chunking was performed in step 1. The ligation step is computationally light and uses variants in the intersection of the chunks to provide chromosome-wide haplotypes. The result of this step (or the previous step if no chunking was used), is used as a haplotype scaffold for the next step.

3. Phase rare variants (MAF < 0.1%) of a chromosome using **phase_rare**. To do this, we use the haplotype scaffold generated in step 2 (or 1) and we proceed in relatively small chunks (e.g. 5Mb) to run many small jobs in parallel. At the end of this step, we have several fully phased chunks across the chromosome.

4. Concatenate the phased chunks generated in step 3 using **bcftools concat -n**. As in the previous step haplotypes have been phased onto a haplotype scaffold, there is no need to ligate the chunks, and the files can be concatenated without decompression and recompression. This makes this step almost instantaneous, even for large cohorts.

This pipeline should be applied on large sample sizes, usually exceeding N=2,000 samples. For smaller sample sizes, just run phase_common. 

## Phasing simulated data
For this tutorial, we work on a dense 10Mb region simulated with msprime for 20,000 European samples. The data is perfectly phased, therefore we will be able to perform haplotype phasing on a large cohort and also to verify the accuracy of the phasing process.

### STEP1: Phasing common variants
SHAPEIT5 phases common variants using the **phase_common** tool, which has been built upon SHAPEIT4. As an input, phase_common requires an unphased file (with AC and AN tags filled up), and automatically sub-sets the file to the desired MAF (e.g. 0.001).

There are different strategies to phase common variants. The first, is to phase the whole chromosome in a single job. This is feasible for SNP array data, but it is usually not optimal for WGS data. Therefore we recommend to chunk the chromosome into large chunks (e.g. 20 cM) if using large WGS data. Of note, we provinde 20cM chunks for the b38 built in the resource folder. The chunking used in this tutorial is under optimal for real data, such as the UK Biobank dataset. For that please refer to the UKB tutorials for more advanced examples.

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


The simulated data that we use in this tutorial are located in `test/wgs`.

#### Option1: Phasing all common variants in one chunk
SHAPEIT5 can phase common variants (MAF >= 0.001) for the the whole 10Mb region using the following command:

<div class="code-example" markdown="1">
```bash
phase_common --input wgs/target.unrelated.bcf --filter-maf 0.001 --region 1 --map info/chr1.gmap.gz --output tmp/target.scaffold.bcf --thread 8
```
</div>

As this is run on the entire region of the simulation, **no ligation step is required** here. On real data, a ligation step will usually be necessary.

#### Option2: Phasing common variants in multiple chunks and ligate
More realistically, when using WGS data on large sample size, it is good practice to run **phase_common** in different large regions of the chromosomes (e.g. 20cM). In the following, we perform phasing in two overlapping 6Mb windows to showcase ligation. The intersection of the two chunks is therefore 2Mb, that is large enough to have a good amount of heterozygous sites for the ligation step.

<div class="code-example" markdown="1">
```bash
#first chunk: 1:1-6000000
phase_common --input wgs/target.unrelated.bcf --filter-maf 0.001 --region 1:1-6000000 --map info/chr1.gmap.gz --output tmp/target.scaffold.chunk0.bcf --thread 8
phase_common --input wgs/target.unrelated.bcf --filter-maf 0.001 --region 1:4000001-10000000 --map info/chr1.gmap.gz --output tmp/target.scaffold.chunk1.bcf --thread 8
```
</div>

This approach allows for parallel runs of different chunks of data. However, as we phasing of rare variants needs a whole-chromosome haplotype scaffold, we need to ligate these chunks to create a single file.

Ligation of multiple chunks of haplotypes is performed using the **ligate** program. The program requires an ordered list of chunks. We recommend using appropriate naming of the files in the previous step, so that a command such as `ls -1v` can directly produce the list of files in the right order.

<div class="code-example" markdown="1">
```bash
ls -1v tmp/target.scaffold.chunk*.bcf > tmp/files.txt
ligate --input tmp/files.txt --output tmp/target.scaffold.bcf --thread 8 --index
```
</div>

At the end of this step, we have a chromosome-wide haplotype scaffold in the file `tmp/target.scaffold.bcf`.

### STEP2: Phase rare variants small region
We can now proceed with the second stage of the process: phasing rare variants iteratively on the scaffold. This is accomplished by the **phase_rare** program. This needs to be done in chunks, in order to maximise the parallelization of our jobs. In this example, we use ~2.5Mb chunks of data with 0.5Mb of scaffold data on each side. For simplicity, we provide chunks coordinates in the file `info/chunks.coordinates.txt`. The program requires two input files, one is the haplotype scaffold generated in the previous steps, and other input require is a file containing the whole unphased region, the same file provided to **phase_common**: SHAPEIT5 automatically discards common variants for this file, duplicated in the phased scaffold. To process the whole dataset, it's enough to loop over chunk and phase them each one at a time.

<div class="code-example" markdown="1">
```bash
while read LINE; do
	ID=$(echo $LINE | awk '{ print $1; }')
	SRG=$(echo $LINE | awk '{ print $3; }')
	IRG=$(echo $LINE | awk '{ print $4; }')
	phase_rare --input wgs/target.unrelated.bcf --scaffold tmp/target.scaffold.bcf --map info/chr1.gmap.gz --input-region $IRG --scaffold-region $SRG --output tmp/target.phased.chunk$CHK\.bcf  --thread 8
done < info/chunks.coordinates.txt
```
</div>

All these jobs can be run in parallel. By default, **phase_rare** discards buffer regions and only outputs chunks that do not overlap with others. Therefore, at the end of this step, we can simply concatenate all the chunks together in a straightforward way.

### STEP3: Obtaining chromosome-wide phased data
Here, we can concatenate  all files without decompression/recompression. This is performed very quickly with [bcftools concat \-\-naive](https://samtools.github.io/bcftools/bcftools.html#concat). We again recommend naming the files in an appropriate way, so that a command such as `ls -1v` can directly produce the list of files in the right order.

<div class="code-example" markdown="1">
```bash
ls -1v tmp/target.phased.chunk$CHK\.bcf > tmp/files.txt
bcftools concat --naive -f tmp/files.txt -o target.phased.bcf --threads 8
bcftools index -f target.phased.bcf
```
</div>
