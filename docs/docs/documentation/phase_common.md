---
layout: default
title: phase_common
nav_order: 1
parent: Documentation
---
# phase_common
{: .no_toc .text-center }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

### Description
Tool to phase common sites, typically SNP array data, or the first step of WES/WGS data.
This is basically an improved version of [SHAPEIT4](https://odelaneau.github.io/shapeit4/) and serves as a replacement.
If you want to prephase your SNP array data for a GWAS, this is the tool yu're supposed to use.

### Usage1: simple phasing run of common variants in simulated data
Phasing simulated genotype data in test/10k.

<div class="code-example" markdown="1">
```bash
phase_common --input test/10k/msprime.nodup.bcf --filter-maf 0.001  --output test/10k/msprime.common.phased.bcf --region 1 --thread 8
```
</div>
The program phases common variants (\-\-filter-maf 0.001) from the input file (\-\-input test/10k/msprime.nodup.bcf) using 8 threads (\-\-thread 8) on the full chromosome 1 (\-\-region 1) and saves the results in the output file (\-\-output test/10k/msprime.common.phased.bcf).
No genetic map is used, a recombination rate of 1 cM/Mb is assumed by default. 

---

### Usage2: simple phasing run of real data
Phasing real genotype data in test/1000G.

<div class="code-example" markdown="1">
```bash
phase_common --input test/1000G/unphased.bcf --map test/1000G/chr20.b37.gmap.gz --output test/1000G/phased.bcf --region 20 --thread 8
```
</div>

The program phases all variants from the input file (\-\-input test/1000G/unphased.bcf) using 8 threads (\-\-thread 8) on the full chromosome 20 (\-\-region 20) using hapmap genetic map (\-\-map test/1000G/chr20.b37.gmap.gz) and saves the results in the output file (\-\-output test/1000G/phased.bcf).
Of note, HapMap genetic maps for builds b37 and b38 of the genome can be found in the resources folder. 

---

### Usage3: advanced phasing run of real data
Phasing real genotype data in test/1000G using a reference panel and a scaffold.

<div class="code-example" markdown="1">
```bash
phase_common --input test/1000G/unphased.bcf --map test/1000G/chr20.b37.gmap.gz --reference test/1000G/reference.bcf --scaffold test/1000G/scaffold.bcf --output test/1000G/phased.bcf --region 20 --thread 8
```
</div>

Same as run 2, using also a reference panel of haplotypes (\-\-reference test/1000G/reference.bcf) and a haplotype scaffold for some of the target samples (\-\-scaffold test/1000G/scaffold.bcf). 

---

### Usage4: advanced phasing run of real data using pedigree information
Phasing real genotype data in test/1000G using pedigree information scaffold.

<div class="code-example" markdown="1">
```bash
phase_common --input test/1000G/pedigree.bcf --pedigree test/1000G/pedigree.fam --map test/1000G/chr20.b37.gmap.gz --output test/1000G/phased.bcf --region 20 --thread 8
```
</div>

Same as run 2, using pedigree information for data in input (\-\-pedigree test/1000G/pedigree.fam). 
Kids will be automatically scaffolded at all sites that can be resolved using Mendel logic.
Parents will be statistically phased without any constraints.  

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


#### MCMC parameters

| Option name 	      | Argument| Default              | Description |
|:--------------------|:--------|:---------------------|:-------------------------------------|
| \-\-mcmc-iterations | STRING  | 5b,1p,1b,1p,1b,1p,5m | Iteration scheme of the MCMC (burnin=b, pruning=p, main=m) |
| \-\-mcmc-prune      | FLOAT   | 0.999                | Pruning threshold for genotype graphs (internal memory structures)  |
| \-\-mcmc-noinit     | NA      | NA                   | If specified, phasing initialization by PBWT sweep is disabled |

#### PBWT parameters

| Option name 	      | Argument|  Default  | Description |
|:--------------------|:--------|:----------|:-------------------------------------|
| \-\-pbwt-modulo     | FLOAT   | 0.1       | Storage frequency of PBWT indexes in cM |
| \-\-pbwt-depth      | INT     | 4         | Depth of PBWT indexes to condition on  |
| \-\-pbwt-mac        | INT     | 5         | Minimal Minor Allele Count at which PBWT is evaluated |
| \-\-pbwt-mdr        | FLOAT   | 0.1       | Maximal Missing Data Rate at which PBWT is evaluated |
| \-\-pbwt-window     | INT     | 4         | Run PBWT selection in windows of this size |

#### HMM parameters

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
