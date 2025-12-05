---
layout: default
title: phase_common [SHAPEIT4]
nav_order: 1
parent: Documentation
---
# phase_common [formerly known as SHAPEIT4]
{: .no_toc .text-center }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

### Description
Tool to phase common variant sites, typically SNP array data (i.e. prephasing).
This tool is also used as the first step for phasing large WES/WGS datasets.
This is an improved version of [SHAPEIT4](https://odelaneau.github.io/shapeit4/) and should serve as a replacement.

### Usage1: phasing unrelated samples
Go in the `test` folder and run:

<div class="code-example" markdown="1">
```bash
phase_common --input array/target.unrelated.bcf --region 1 --map info/chr1.gmap.gz --output target.phased.bcf --thread 8
```
</div>

The program phases data in the input (unphased) file (\-\-input `array/target.unrelated.bcf`):
- using 8 threads (\-\-thread 8)
- data located on chromosome 1 (\-\-region 1) 
- using a specific genetic map (\-\-input `info/chr1.gmap.gz`)
- saving the phased haplotypes in output file (\-\-output `target.phased.bcf`).

If no genetic map is specified, a recombination rate of 1 cM/Mb is assumed by default.

---

### Usage2: phasing related samples
Go in the `test` folder and run:

<div class="code-example" markdown="1">
```bash
phase_common --input array/target.family.bcf --pedigree info/target.family.fam --region 1 --map info/chr1.gmap.gz --output target.phased.bcf --thread 8
```
</div>

The family relationships (parent-offspring relationships, duos/trios) are specified in the file `info/target.family.fam`.
This file contains one line per sample having parent(s) in the dataset and three columns (kidID fatherID and motherID), separated by TABs for spaces.
Use NAs for unknown parents (in the case of duos). In output file, the first offspring haplotype is transmitted by the father, the second by the mother.

This tools uses family data to fix the phase of offspring heterozygous genotypes when possible, that is when:
- There is no Mendel inconsistency,
- At least one parent is homozygous.
In other words, it builds a scaffold of haplotypes for offsprings from the parental genomes.   

---

### Usage3: phasing chromosome X data
Go in the `test` folder and run:

<div class="code-example" markdown="1">
```bash
phase_common --input array/target.haploid.bcf --haploids info/target.haploid.txt --region 1 --map info/chr1.gmap.gz --output target.phased.bcf --thread 8
```
</div>

The list of samples being haploid (e.g. males for chromosome X) is specified in the file `info/target.haploid.fam`. This file contains one line per sample.

For convenience, haploid samples are encoded either internally or in input/output as fully homozgote diploid samples.
In practice, this tool sets to missing any heterozygous genotype in input and re-impute them as REF/REF or ALT/ALT.
For all samples not listed in the haploid file, the tool proceeds with standard phasing. 

---

### Usage4: phasing using a reference panel
Go in the `test` folder and run:

<div class="code-example" markdown="1">
```bash
phase_common --input array/target.unrelated.bcf --reference array/reference.bcf --region 1 --map info/chr1.gmap.gz --output target.phased.bcf --thread 8
```
</div>

The tool uses the haplotypes specified in the BCF/VCF file `array/reference.bcf` as a reference panel to phase all the data in input file.
Any input site not in the reference panel is not considered in the phasing.

---

### Usage5: phasing using a scaffold
Go in the `test` folder and run:

<div class="code-example" markdown="1">
```bash
phase_common --input array/target.unrelated.bcf --scaffold array/target.scaffold.bcf --region 1 --map info/chr1.gmap.gz --output target.phased.bcf --thread 8
```
</div>

The tool uses the haplotypes specified in the BCF/VCF file `array/target.scaffold.bcf` as a scaffold to phase all the data in input file.
The scaffold can contain data for a subset of samples at a subset of variant sites. Data in the scaffold needs to be non-missing and phased.

---

### Command line options

#### Basic options

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-\-help             | NA      | NA       | Produces help message |
| \-\-seed             | INT     | 15052011 | Seed of the random number generator  |
| \-T \[ \-\-thread \] | INT     | 1        | Number of thread used|

#### Input files

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-I \[\-\-input \]   | STRING  | NA       | Genotypes to be phased in VCF/BCF/XCF format |
| \-H \[\-\-reference \]| STRING  | NA       | Reference panel of haplotypes in VCF/BCF/XCF format  |
| \-S \[\-\-scaffold \]| STRING  | NA       | Scaffold of haplotypes in VCF/BCF/XCF format  |
| \-M \[\-\-map \]     | STRING  | NA       | Genetic map  |
| \-\-pedigree         | STRING  | NA       | Pedigree information (offspring father mother triplets) |
| \-\-haploids         | STRING  | NA       | List of samples that are haploids (e.g. males for chrX) |
| \-R \[\-\-region \]  | STRING  | NA       | Target region  |


#### Filter parameters

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-\-filter-snp       | NA      | NA       | If specified, the program only consider SNPs |
| \-\-filter-maf       | FLOAT   | 0        | \[Expert option\] Only consider variants with MAF above the specifed value. It requires AC/AN tags in VCF/BCF file. |


#### MCMC parameters [Expert]

| Option name 	      | Argument| Default              | Description |
|:--------------------|:--------|:---------------------|:-------------------------------------|
| \-\-mcmc-iterations | STRING  | 5b,1p,1b,1p,1b,1p,5m | Iteration scheme of the MCMC (burnin=b, pruning=p, main=m) |
| \-\-mcmc-prune      | FLOAT   | 0.999                | Pruning threshold for genotype graphs (internal memory structures)  |
| \-\-mcmc-noinit     | NA      | NA                   | If specified, phasing initialization by PBWT sweep is disabled |

#### PBWT parameters [Expert]

| Option name 	      | Argument|  Default  | Description |
|:--------------------|:--------|:----------|:-------------------------------------|
| \-\-pbwt-modulo     | FLOAT   | 0.1       | Storage frequency of PBWT indexes in cM |
| \-\-pbwt-depth      | INT     | 4         | Depth of PBWT indexes to condition on  |
| \-\-pbwt-mac        | INT     | 5         | Minimal Minor Allele Count at which PBWT is evaluated |
| \-\-pbwt-mdr        | FLOAT   | 0.1       | Maximal Missing Data Rate at which PBWT is evaluated |
| \-\-pbwt-window     | INT     | 4         | Run PBWT selection in windows of this size |

#### HMM parameters [Expert]

| Option name 	      | Argument|  Default  | Description |
|:--------------------|:--------|:----------|:-------------------------------------|
| \-\-hmm-window      | INT     | 4         | Minimal size of the phasing window in cM |
| \-\-hmm-ne          | INT     | 15000     | Effective size of the population |

#### Output files

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-O \[\-\-output \]  | STRING  | NA       | Phased haplotypes in VCF/BCF/XCF format |
| \-\-output-format     | STRING  | bcf       | File format for the output ([bcf] standard VCF/BCF format / [graph] graph format that intergrates phasing uncertainty / [bh] XCF binary format for fast loading in Impute5)  |
| \-\-log              | STRING  | NA       | Log file  |
