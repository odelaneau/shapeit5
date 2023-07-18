---
layout: default
title: UK Biobank SNP array data
nav_order: 2
parent: Tutorials
---
# UK Biobank SNP array data
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

**For any question on this pipeline, please contact [Robin J. Hofmeister](https://github.com/RJHFMSTR) and [Olivier Delaneau](https://github.com/odelaneau).** 

---

**IMPORTANT: For this tutorial, we work on the UK Biobank Research Analysis Platform (RAP).**

## Rationale
Although designed to phase rare variants from sequencing data, SHAPEIT5 contains several tools and therefore can also be used to phase common variant only, for example from SNP array technologies. This can be done using a single tool: **SHAPEIT5_phase_common**. As the number of marker on a SNP array is much small than in WGS data, the phasing of SNP array can be done in a single job per chromosome.

In this tutorial, we cover the phasing of the UK Biobank Axiom array across ~500k individuals on the UKB RAP platform.




<br>
## Phasing of SNP array data
<br>
### Set up your environment
To structure the outputs of the analysis and to simplify the understanding of each step of this tutorial, we first create output folders as follows. You can choose the change the name of these folders but you will have to change the code accordingly.
<div class="code-example" markdown="1">
```bash
dx mkdir -p Phasing/PhasingSNParray/step1_dataqc/
dx mkdir -p Phasing/PhasingSNParray/step2_chrrename/
dx mkdir -p Phasing/PhasingSNParray/step3_swapalleles/
dx mkdir -p Phasing/PhasingSNParray/step4_liftover/
dx mkdir -p Phasing/PhasingSNParray/step5_phasing/
```
</div>
<br>
### Quality Control
For the phasing of the UK Biobank SNP array data, we first perform a quality control of the data. For this, we use the UK Biobank SNPs and samples QC file (UK Biobank Resource 531) to only retain SNPs and individuals that have been used for the official phasing of the Axiom array data.

For your information, this file can be found at this location on the UKB RAP:  `/Bulk/Genotype Results/Genotype calls/ukb_snp_qc.txt`. You don't need to download it, it is automatically located within our scripts.

To quality control the SNP array data, proceed as follows, which will produce the output file `Phasing/PhasingSNParray/step1_dataqc/full_c${CHR}_b0_v2.vcf.gz`
<div class="code-example" markdown="1">
```bash
#Get SNP list
dx run app-swiss-army-knife -iin="/Bulk/Genotype Results/Genotype calls/ukb_snp_qc.txt" --folder="/Phasing/PhasingSNParray/step1_dataqc/" -icmd="cat ukb_snp_qc.txt | awk '{ print \$1, \$159; }' > SNPlist.unfiltered.txt && cat SNPlist.unfiltered.txt | sed '1d' | awk '{ if (\$2 == 1) print \$1; }' > SNPlist.filtered.QC.txt" --instance-type mem1_ssd1_v2_x2 --name qc_snp --priority normal -y

#Get sample lists
dx run app-swiss-army-knife --folder="/Phasing/PhasingSNParray/step1_dataqc/" -iin="Bulk/Genotype\ Results/Genotype\ calls/ukb_sqc_v2.txt" -iin="Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c1_b0_v2.fam" -icmd="cat ukb_sqc_v2.txt | cut -d ' ' -f 66 > in_phasing.txt && cat ukb22418_c1_b0_v2.fam | cut -d ' ' -f 1 > samples.tmp.txt && paste samples.tmp.txt in_phasing.txt > INDlist.unfiltered.txt && cat INDlist.unfiltered.txt | sed '1d' | awk '{ if (\$2 == 1) print \$1, \$1; }' > INDlist.filtered.QC.txt && rm samples.tmp.txt && rm in_phasing.txt" --instance-type mem1_ssd1_v2_x4 --priority normal --name qc_sample -y


#Filter each chromosome
for CHR in $(seq 1 22); do
	dx run app-swiss-army-knife -iin="/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c${CHR}_b0_v2.bed" -iin="/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c${CHR}_b0_v2.bim" -iin="/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c${CHR}_b0_v2.fam" -iin="/Phasing/PhasingSNParray/step1_dataqc/INDlist.filtered.QC.txt" -iin="/Phasing/PhasingSNParray/step1_dataqc/SNPlist.filtered.QC.txt" --folder="/Phasing/PhasingSNParray/step1_dataqc/" -icmd="plink2 --bfile ukb22418_c${CHR}_b0_v2 --keep INDlist.filtered.QC.txt --extract SNPlist.filtered.QC.txt --export vcf bgz --out full_c${CHR}_b0_v2 && bcftools index full_c${CHR}_b0_v2.vcf.gz" --instance-type mem1_ssd1_v2_x2 --name qc_chr${CHR} --priority normal -y
done
```
</div>

<br>


### Lift over (optional)
To match the UKB SNP array data (hg19) to the UKB sequencing data (hg38), we lift over the UKB SNP array data from hg19 to hg38. This step is not necessary if you only want to phase the SNP array data. In our case, this allow us to merge the SNP array data with the UKB whole-exome sequencing data (WES) to phase the WES data, as described in the dedicated tutorial.

**Requirement:**
This step relies on two additional codes (`swaprefalt_0.0.1.tar.gz` and `liftovervcf_0.0.1.tar.gz`) that are provided as docker images [here](https://github.com/odelaneau/shapeit5/releases). Make sure to upload these docker images in your `/docker` folder on the UKB RAP using:

<br>
<div class="code-example" markdown="1">
```bash
# download required docker images from our github repository
wget https://github.com/odelaneau/shapeit5/releases/download/v1.0.0/swaprefalt_0.0.1.tar.gz
wget https://github.com/odelaneau/shapeit5/releases/download/v1.0.0/liftovervcf_0.0.1.tar.gz

# upload the docker images on the UKB RAP
dx mkdir -p docker/
dx upload swaprefalt_0.0.1.tar.gz --path="docker/"
dx upload liftovervcf_0.0.1.tar.gz --path="docker/"
```
</div>
<br>

The lift over we perform consist of three phases:

<div class="code-example" markdown="1">
- Rename chromosomes
- Swap alleles
- Lift Over
</div>


<br>

#### Rename chromosomes

<div class="code-example" markdown="1">
```bash
# create a file to match chromosome tag between hg19 and hg38 version
for CHR in {1..22}; do echo "${CHR} chr${CHR}"; done > chr_rename.txt

# upload this file on the UKB RAP
dx upload chr_rename.txt --path="Phasing/PhasingSNParray/step2_chrrename/"

# rename chromosomes
ANN=/mnt/project/Phasing/PhasingSNParray/step2_chrrename/chr_rename.txt
for CHR in {1..22}; do
	VCFF=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/full_c$CHR\_b0_v2.vcf.gz
	OUTF=full_c$CHR\_b0_v2.b37.vcf.gz
	dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step2_chrrename/" -icmd="bcftools annotate -Oz -o $OUTF --rename-chrs $ANN $VCFF && bcftools index $OUTF" --instance-type mem2_ssd1_v2_x2 --name updatechr_chr$CHR --priority normal -y
	
done
```
</div>




#### Swap alleles

<div class="code-example" markdown="1">
```bash
for CHR in 20; do
	VCFF=/mnt/project/Phasing/PhasingSNParray/step2_chrrename/full_c$CHR\_b0_v2.b37.vcf.gz
	OUTF=full_c$CHR\_b0_v2.b37.swapped.vcf.gz
	dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step3_swapalleles/" -iimage_file="/docker/swaprefalt_0.0.1.tar.gz" -icmd="swapRefAlt_static --input $VCFF --output $OUTF && bcftools index $OUTF" --instance-type mem2_ssd1_v2_x2 --priority normal --name swapalleles_chr$CHR -y
done
```
</div>

#### Lift over

<div class="code-example" markdown="1">
```bash
# Download the reference fasta
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
dx upload GRCh38_full_analysis_set_plus_decoy_hla.fa --path="Phasing/PhasingSNParray/step4_liftover/"

# Download the chain file
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
dx upload hg19ToHg38.over.chain.gz --path="Phasing/PhasingSNParray/step4_liftover/"

for CHR in {1..22}; do
	CHAIN=/mnt/project/Phasing/PhasingSNParray/step4_liftover/hg19ToHg38.over.chain.gz
	VCFF=/mnt/project/Phasing/PhasingSNParray/step3_swapalleles/full_c$CHR\_b0_v2.b37.swapped.vcf.gz
	REF=/mnt/project/Phasing/PhasingSNParray/step4_liftover/GRCh38_full_analysis_set_plus_decoy_hla.fa
	LIFF=full_c$CHR\_b0_v2.b38.vcf.gz
	SORF=full_c$CHR\_b0_v2.b38.sorted.vcf.gz
	dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step4_liftover/" -iimage_file="/docker/liftovervcf_0.0.1.tar.gz" -icmd="liftoverVCF_static --input $VCFF --output $LIFF --chain $CHAIN --fasta $REF --chr chr$CHR && bcftools sort -Oz -m 6G -o $SORF $LIFF && rm $LIFF && bcftools index $SORF" --instance-type mem2_ssd1_v2_x2 --priority normal --name liftover_chr$CHR -y
done
```
</div>



<br>

### Phasing 

**IMPORTANT**: in the following code make sure to change the shapeit5 docker image name (here `shapeit5_beta.tar.gz`) to the latest version that you've downloaded [here](https://odelaneau.github.io/shapeit5/docs/installation/docker)
.
<div class="code-example" markdown="1">
```bash
# Download map files
dx mkdir -p data/shapeit_maps/
wget https://github.com/odelaneau/shapeit5/raw/main/maps/genetic_maps.b38.tar.gz
tar -xvzf genetic_maps.b38.tar.gz
dx upload *.b38.gmap.gz --path="data/shapeit_maps/"

# Phasing
for CHR in {1..22}; do
	MAP=/mnt/project/data/shapeit_maps/chr${CHR}.b38.gmap.gz
	INP=/mnt/project/Phasing/PhasingSNParray/step4_liftover/full_c$CHR\_b0_v2.b38.sorted.vcf.gz
	LOG=c${CHR}_b0_v2.b38.sorted.phased.log
	OUT=c${CHR}_b0_v2.b38.sorted.phased.bcf
	TIM=c${CHR}_b0_v2.b38.sorted.phased.time

	dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step5_phasing/" -iimage_file="/docker/shapeit5_beta.tar.gz" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase_common --input $INP --map $MAP --output $OUT --region chr${CHR} --log $LOG --thread 16 && bcftools index $OUT --threads 16" --instance-type mem2_ssd1_v2_x16 --priority normal --name phasing_chr${CHR} -y
done

```
</div>


The full list of options for the **SHAPEIT5_phase_common** command can be found [here](https://odelaneau.github.io/shapeit5/docs/documentation/phase_common/).


<br>
## Validation of your phasing
You can validate the quality of the haplotypes using the **SHAPEIT5_switch** tool. For this you will need parent-offspring duos or trios stored in a three-columns file (here called `family.ped`) following this format:

<div class="code-example" markdown="1">
- family.ped:      `offspring_id`    `parent1_id`    `parent2_id`

</div>


To validate the phasing using family data, the phasing must be performed by excluding parental genomes, so that offsprings are phased regardless of their parental genomes. This can be done using the **bcftools view** command, with as input the original SNP array data after the quality control and liftover steps (located here `/Phasing/PhasingSNParray/step4_liftover/full_c$CHR\_b0_v2.b38.sorted.vcf.gz`)

<div class="code-example" markdown="1">
```bash
dx mkdir -p Phasing/PhasingSNParray/benchmark/
for CHR in {1..22}; do
IN=/mnt/project/Phasing/PhasingSNParray/step4_liftover/full_c$CHR\_b0_v2.b38.sorted.vcf.gz
OUT=benchmark_c$CHR\_b0_v2.b38.sorted.vcf.gz
dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/benchmark/" -icmd="bcftools view --threads 16 -S ^parents.txt -Ob -o ${OUT} ${IN} && bcftools index ${OUT} --threads 16" --instance-type mem1_ssd1_v2_x16 --priority normal --name phasing_chr${CHR} -y
done
```
</div>

After excluding parental genomes with the above command, proceed with the normal phasing procedure as detailed above.

Let's consider that you performed the above steps of phasing using the input data exluding parental genomes, which produced a phased output file that you named `benchmark_c${CHR}_b0_v2.b38.sorted.phased.bcf`). You can validate your phasing using the following command:


<div class="code-example" markdown="1">
```bash
for CHR in {1..22}; do
dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/benchmark/" -iimage_file="/docker/shapeit5_beta.tar.gz" -icmd="SHAPEIT5_switch --validation /Phasing/PhasingSNParray/step4_liftover/full_c$CHR\_b0_v2.b38.sorted.vcf.gz --estimation benchmark_c${CHR}_b0_v2.b38.sorted.phased.bcf --region ${CHR} --output benchmark_chr${CHR}"  --instance-type mem2_ssd1_v2_x16 --priority normal --name benchmark_chr${CHR} -y
done
```
</div>



The full list of options for the **SHAPEIT5_switch** command can be found [here](https://odelaneau.github.io/shapeit5/docs/documentation/switch/).





















