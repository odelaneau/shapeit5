---
title: Expert quick start
nav_order: 10
---

# Quick Start  [for experts]

If you are familiar with phasing and the pipeline of SHAPEIT5, appling SHAPEIT5 on your data will be pretty straightforward to understand.
Instructions below provide quick steps for phasing your data.


## Phase_common
This command phase only common variants (above the specified MAF threshold). The output of this phasing is used as a scaffold to phase rare variants. If this phasing is done by chunks, you will have to ligate the chunks using **SHAPEIT5_ligate* tool. If this phasing is done for entire chromosome, you don't need to ligate.

<div class="code-example" markdown="1">
```bash
SHAPEIT5_phase_common --input $IN --map $MAP --output $OUT --thread $THREADS --log $LOG --filter-maf 0.001 --region $REG
```
</div>

## Ligate
This command assemble the phasing chunks to keep the correct phase while assembling. Typically, it take as input the list of chunks produces with **SHAPEIT5_phase_common**.
<div class="code-example" markdown="1">
```bash	
SHAPEIT5_ligate --input $IN --output $OUT --thread $THREADS --index
```
</div>



## Phase_rare
This command phase rare variants onto the haplotype scaffold of common variants. Typically, the scaffold is the ligated chunks from **SHAPEIT5_ligate**.

<div class="code-example" markdown="1">
```bash	
SHAPEIT5_phase_rare --input-plain $IN --map $MAP --scaffold $SCAF --output $OUT --thread $THREADS --log $LOG --filter-maf 0.001 --region $REG
```
</div>


