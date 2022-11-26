---
layout: default
title: phase_rare
nav_order: 3
parent: Documentation
---
# phase_rare
{: .no_toc .text-center }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

### Description
Tool to phase rare variants onto a scaffold of common variants (output of phase_common / ligate).

### Usage
Simple run

<div class="code-example" markdown="1">
```bash
#chunk1
SHAPEIT5_phase_rare --input-plain 10k/msprime.nodup.bcf --scaffold test/10k/msprime.common.truth.bcf --output test/10k/msprime.rare.chunk1.bcf --scaffold-region 1:1000000-3000000 --input-region 1:1500000-2500000 --thread 8

#chunk2
SHAPEIT5_phase_rare --input-plain 10k/msprime.nodup.bcf --scaffold test/10k/msprime.common.truth.bcf --output test/10k/msprime.rare.chunk2.bcf --scaffold-region 1:2000000-4000000 --input-region 1:2500001-3500000 --thread 8
```
</div>

The first command phases rare variants from the input file (\-\-input-plain test/10k/msprime.nodup.bcf) using 8 threads (\-\-thread 8) on the region 1500000-2500000 of chromosome 1 (\-\-input-region 1:1500000-2500000) using a haplotype scaffold phased for the full chromosome 1 obtained from the phase_common program in the region 1000000-3000000 of chromosome 1 (\-\-scaffold-region 1:1000000-3000000) and saves the results in the output file (\-\-output test/10k/msprime.rare.chunk1.bcf).

The obtained files can be quickly concatenated to generate chromosome-wide files using bcftools concat --naive.

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
| \-\-input-plain      | STRING  | NA       | Genotypes to be phased in plain VCF/BCF format |
| \-\-input-region     | STRING  | NA       | Region to be considered in \-\-input-plain |
| \-\-input-maf        | FLOAT   | 0.001    | Threshold for sparse genotype representation in --input-plain |
| \-\-scaffold         | STRING  | NA       | Scaffold of haplotypes in VCF/BCF format  |
| \-\-scaffold-region  | STRING  | NA       | Region to be considered in \-\-scaffold  |
| \-\-map              | STRING  | NA       | Genetic map  |
| \-\-pedigree         | STRING  | NA       | Pedigree information (offspring father mother) |


#### PBWT parameters

| Option name 	      | Argument|  Default  | Description |
|:--------------------|:--------|:----------|:-------------------------------------|
| \-\-pbwt-modulo     | FLOAT   | 0.1       | Storage frequency of PBWT indexes in cM |
| \-\-pbwt-depth-common | INT     | 2         | Depth of PBWT indexes at common sites to condition on  |
| \-\-pbwt-depth-rare | INT     | 2         | Depth of PBWT indexes at rare sites to condition on  |
| \-\-pbwt-mac        | INT     | 2         | Minimal Minor Allele Count at which PBWT is evaluated |
| \-\-pbwt-mdr        | FLOAT   | 0.1       | Maximal Missing Data Rate at which PBWT is evaluated |

#### HMM parameters

| Option name 	      | Argument|  Default  | Description |
|:--------------------|:--------|:----------|:-------------------------------------|
| \-\-effective-size  | INT     | 15000     | Effective size of the population |

#### Output files

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-O \[\-\-output \]  | STRING  | NA       | Phased haplotypes in VCF/BCF format |
| \-\-output-buffer    | STRING  | NA       | If specified, right and left buffers are printed in output |
| \-\-log              | STRING  | NA       | Log file  |
