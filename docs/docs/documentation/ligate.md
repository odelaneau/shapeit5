---
layout: default
title: ligate
nav_order: 2
parent: Documentation
---
# ligate
{: .no_toc .text-center }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

### Description
Ligate multiple phased BCF/VCF files into a single whole chromosome file. Typically run to ligate multiple chunks of phased common variants.

### Usage
Simple run

<div class="code-example" markdown="1">
```bash
#ls -1v in order to keep the order within the chromosome
ls -1v chr1/*.phased.bcf > list_phased_files_chr1.txt

SHAPEIT5_ligate --input list_phased_files_chr1.txt --output ligated_chr1.bcf --thread 2
```
</div>

The program ligates together multiple phased files, listed in the input txt file (\-\-input list_phased_files.txt) overlapping at buffer regions, using two threads (\-\-thread 2), and saves the ouput file (\-\-output ligated_chr1.bcf).

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
| \-I \[\-\-input \]   | STRING  | NA       | Text file containing all VCF/BCF to ligate, one file per line |

#### Output files

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-O \[\-\-output \]  | STRING  | NA       | Phased haplotypes in VCF/BCF format |
| \-\-no-index         | STRING  | NA       | If specified, the output the the ligated VCF/BCF is not indexed for random access to genomic regions |
| \-\-log              | STRING  | NA       | Log file  |
