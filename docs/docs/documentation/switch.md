---
layout: default
title: switch
nav_order: 4
parent: Documentation
---
# switch
{: .no_toc .text-center }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

### Description
Program to compute switch error rate (SER) and genotyping error rate (GER) given validation data.
Validation data, in which haplotypes are known, can be obtained in multiple ways:
- Simulated data, using msprime for instance.
- Trio/duo data. Estimate haplotypes for offspring excluding the parents from the dataset. Use the parent and the offsprings in the validation set.
- Chromosome X data. Pair male chromosomes X to make females in which phase in known.

### Usage1: Benchmark haplotypes using simulated data

First, run a phasing run:
 
<div class="code-example" markdown="1">
```bash
phase_common --input array/target.unrelated.bcf --region 1 --map info/chr1.gmap.gz --output target.phased.bcf --thread 8
```
</div>

The data in `array/target.unrelated.bcf` is phased simulated data, so we can also use it as validation set.
Open the BCF file to check:
<div class="code-example" markdown="1">
```bash
bcftools view -H array/target.unrelated.bcf | cut -f1-20 | head
```
</div>

Then, run thw switch program to compare the true haplotypes to those estimated by phase_common:
 
<div class="code-example" markdown="1">
```bash
switch --validation array/target.unrelated.bcf --estimation target.phased.bcf --region 1 --output target.phased --thread 8
```
</div>

This command will produce a bunch of files:
- `target.phased.block.switch.txt.gz`has 2 columns (sample ID, position). This gives the coordinates of blocks of data coorectly phased. Can be used to produced a visual representation of the phasing   . See supplemetary figure 1 of [SHAPEIT4 paper](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-13225-y/MediaObjects/41467_2019_13225_MOESM1_ESM.pdf).
- `target.phased.frequency.switch.txt.gz` has 4 columns (MAC, #errors, #hets, SER). This gives SER stratified by MAC.
- `target.phased.sample.switch.txt.gz` has 4 columns (sample ID, #errors, #hets, SER). This gives the SER per sample. This is the most important information given by the switch program.
- `target.phased.sample.typing.txt.gz` has 4 columns (sample ID, #errors, #genotypes, GER). This gives GER per sample.
- `target.phased.variant.switch.txt.gz` has 5 columns (rsid, position, #errors, #hets, SER). This gives the SER per variant relative to previous hets. 
- `target.phased.variant.typing.txt.gz` has 5 columns (rsid, position, #errors, #genotypes, GER). This gives the GER per variant. 

To compute the SER across the entire dataset, we recommend to sum the numbers of errors and hets across all samples first:
<div class="code-example" markdown="1">
```bash
zcat target.phased.sample.switch.txt.gz | awk 'BEGIN { e=0; t=0; } { e+=$2; t+=$3; } END { print "SER =", e*100/t; }'
```
</div>

### Usage2: Benchmark haplotypes using family data

First, run build a benchmark dataset by removing parental genomes:
<div class="code-example" markdown="1">
```bash
cat info/target.family.fam | cut -f2- | tr "\t" "\n" > parents.txt
bcftools view -Ob -o benchmark.data.bcf -S ^parents.txt array/target.family.bcf
bcftools index benchmark.data.bcf
```
</div>

Second, phase the benchmark dataset:
<div class="code-example" markdown="1">
```bash
phase_common --input benchmark.data.bcf --region 1 --map info/chr1.gmap.gz --output target.phased.bcf --thread 8
```
</div>

Third, validate it using family data (original BCF + FAM file):
<div class="code-example" markdown="1">
```bash
switch --validation array/target.family.bcf --estimation target.phased.bcf --pedigree info/target.family.fam --region 1 --output target.phased --thread 8
```
</div>

Fourth, compute SER:
<div class="code-example" markdown="1">
```bash
zcat target.phased.sample.switch.txt.gz | awk 'BEGIN { e=0; t=0; } { e+=$2; t+=$3; } END { print "SER =", e*100/t; }'
```
</div>

### Usage3: Benchmark haplotypes using haploid chromosome X data

To come!

---

### Command line options

#### Basic options

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-\-help             | NA      | NA       | Produces help message |
| \-T \[ \-\-thread \] | INT     | 1        | Number of thread used|

#### Input files

| Option name 	          | Argument| Default  | Description |
|:------------------------|:--------|:---------|:-------------------------------------|
| \-V \[\-\-validation \] | STRING  | NA       | Validation dataset in VCF/BCF format |
| \-E \[\-\-estimation \] | STRING  | NA       | Phased dataset in VCF/BCF format  |
| \-F \[\-\-frequency \]  | STRING  | NA       | Variant frequency in VCF/BCF format, to exaclude variants and/or stratify SER by MAC |
| \-P \[\-\-pedigree \]   | STRING  | NA       | Pedigree information (offspring father mother) |
| \-R \[\-\-region \]     | STRING  | NA       | Target region  |
| \-\-nbins               | INT     | 20       | Number of bins used for calibration (for PP field) |
| \-\-min-pp              | FLOAT   | 0        | Minimal PP value for entering computations |
| \-\-singleton           | STRING  | NA       | Singleton phase |
| \-\-dupid               | STRING  | NA       | Duplicate ID for UKB matching IDs |

#### Output files

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-O \[\-\-output \]  | STRING  | NA       | Phased haplotypes in VCF/BCF format |
| \-\-log              | STRING  | NA       | Log file  |

