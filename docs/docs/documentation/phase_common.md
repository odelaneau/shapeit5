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

### Usage
Simple run

<div class="code-example" markdown="1">
```bash
SHAPEIT5_phase_common --input test/10k/msprime.nodup.bcf --filter-maf 0.001  --output test/10k/msprime.common.phased.bcf --region 1 --thread 8
```
</div>

The program phases common variants (\-\-filter-maf 0.001) from the input file (\-\-input test/10k/msprime.nodup.bcf) using 8 threads (\-\-thread 8) on the full chromosome 1 (\-\-region 1) and saves the results in the output file (\-\-output test/10k/msprime.common.phased.bcf).

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
| \-I \[\-\-input \]   | STRING  | NA       | Genotypes to be phased in VCF/BCF format |
| \-H \[\-\-reference \]| STRING  | NA       | Reference panel of haplotypes in VCF/BCF format  |
| \-S \[\-\-scaffold \]| STRING  | NA       | Scaffold of haplotypes in VCF/BCF format  |
| \-M \[\-\-map \]     | STRING  | NA       | Genetic map  |
| \-\-pedigree         | STRING  | NA       | Pedigree information (offspring father mother) |
| \-R \[\-\-region \]  | STRING  | NA       | Target region  |


#### Filter parameters

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-\-filter-snp       | NA      | NA       | If specified, the program only consider SNPs |
| \-\-filter-maf       | FLOAT   | 0        | \[Expert option\] Only consider variants with MAF above the specifeed value. It requires AC/AN tags in VCF/BCF file. |


#### MCMC parameters

| Option name 	      | Argument| Default              | Description |
|:--------------------|:--------|:---------------------|:-------------------------------------|
| \-\-mcmc-iterations | STRING  | 5b,1p,1b,1p,1b,1p,5m | Iteration scheme of the MCMC (burnin=b, pruning=p, main=m) |
| \-\-mcmc-prune      | FLOAT   | 0.999                | Pruning threshold for genotype graphs  |
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
| \-O \[\-\-output \]  | STRING  | NA       | Phased haplotypes in VCF/BCF format |
| \-\-output-graph     | STRING  | NA       | Phased haplotypes in BIN format (Useful to sample multiple likely haplotype configurations per sample)  |
| \-\-log              | STRING  | NA       | Log file  |
