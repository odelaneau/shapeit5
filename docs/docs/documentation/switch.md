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
Program to compute switch error rate and genotyping error rate given simulated or trio data.

### Usage
Simple run

<div class="code-example" markdown="1">
```bash
SHAPEIT5_switch --validation 10k/msprime.nodup.bcf --estimation 10k/msprime.rare.chunk1.bcf --region 1 --output 10k/msprime.rare.chunk1
```
</div>

The program estimates errors from the phased file (\-\-estimation 10k/msprime.rare.chunk1.bcf) on the full chromosome 1 (\-\-region 1) using the a validation file \-\-validation 10k/msprime.nodup.bcf) and saves the results in several output files with the specified prefix (\-\-output 10k/msprime.rare.chunk1).

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
| \-F \[\-\-frequency \]  | STRING  | NA       | Variant frequency in VCF/BCF format  |
| \-P \[\-\-pedigree \]   | STRING  | NA       | Pedigree file in PED format  |
| \-R \[\-\-region \]     | STRING  | NA       | Target region  |
| \-\-nbins               | INT     | 20       | Number of bins used for calibration |
| \-\-min-pp              | FLOAT   | 0        | Minimal PP value for entering computations |
| \-\-singleton           | STRING  | NA       | Singleton phase |
| \-\-dupid               | STRING  | NA       | Duplicate ID for UKB matching IDs |

#### Output files

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-O \[\-\-output \]  | STRING  | NA       | Phased haplotypes in VCF/BCF format |
| \-\-log              | STRING  | NA       | Log file  |
