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
Ligate multiple phased BCF/VCF files into a single whole chromosome file. Typically used to ligate multiple chunks of phased common variants.

### Usage: Ligate two chunks of phased data
Go in the `test` folder and run:

First, let's phase two overlapping chunks of data:
<div class="code-example" markdown="1">
```bash
phase_common --input array/target.unrelated.bcf --region 1:1-5500000 --map info/chr1.gmap.gz --output target.chunk1.bcf --thread 8
phase_common --input array/target.unrelated.bcf --region 1:4500000-10000000 --map info/chr1.gmap.gz --output target.chunk2.bcf --thread 8
```
</div>

Then, we make a text file listing the chunks of data to ligate. *This file need to be sorted!*
<div class="code-example" markdown="1">
```bash
echo target.chunk1.bcf > chunks.txt
echo target.chunk2.bcf >> chunks.txt
```
</div>

Finally, we can proceed with the ligation:
<div class="code-example" markdown="1">
```bash
ligate --input chunks.txt --output target.phased.bcf --thread 2
```
</div>

The ligate program uses two threads (\-\-thread 2) and saves the ligated haplotypes in the output file (\-\-output `target.phased.bcf`).
Two statistics are reported while during ligation:
- Switch rate: the percentage of samples for which haplotypes have been switched,
- Avg phaseQ: the agreement of the phasing in the overlapping region between two successive chunks of data. We usually expect values above 80 or 90.

### Note: Ligate chunks with family information

If some samples have been phased using pedigree information (trios/duos, option \-\-pedigree in phase_common), you must use \-\-pedigree here too, with the same ped file, to make sure that these haplotypes are not switched when ligated.
If this is not done, haplotypes may be switched for offsprings which would make the paternal/maternal labels incorrects.

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
| \-\-pedigree         | STRING  | NA       | Pedigree information (offspring father mother triplets). Really important to make sure that scaffolded samples are not switched! |

#### Output files

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-O \[\-\-output \]  | STRING  | NA       | Phased haplotypes in VCF/BCF format |
| \-\-no-index         | STRING  | NA       | If specified, the output the the ligated VCF/BCF is not indexed for random access to genomic regions |
| \-\-log              | STRING  | NA       | Log file  |
