---
layout: default
title: UK Biobank WGS data
nav_order: 3
parent: Tutorials
---
# UK Biobank WGS data
{: .no_toc }

{: .warning }
Website under construction: content not available yet!

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---


## Rational
SHAPEIT5 is a two-step approach that treats each chromosome independently and works as follows:

1. Phase common variants (MAF >= 0.1%) of a chromosome using SHAPEIT5_phase_common. This can be done as a single job for SNP array or WES data, but it might be necessary to split the chromosome into large chunks (e.g. 20 cM) in the case of WGS data.

2. Ligate the phased common variants (MAF >= 0.1%) of a chromosome using SHAPEIT5_ligate, only if chunking was performed in step 1. The ligation step is computationally light and uses variants in the intersection of the chunks to provide chromosome-wide haplotypes. The result of this step (or the previous step if no chunking was used), is used as a haplotype scaffold for the next step.

3. Phase rare variants (MAF < 0.1%) of a chromosome using SHAPEIT5_rare. To do this, we use the haplotype scaffold generated in step 2 (or 1) and we proceed in relatively small chunks (e.g. 5Mb) to run relatively fast job in parallel. At the end of this step, we have several fully phased chunks across the chromosome.

4. Concatenate the phased chunks generated in step 3 using bcftools concat -n. As in the previous step haplotypes have been phased onto a haplotype scaffold, there is no need to ligate the chunks, and the files can be concatenated without decompression and recompression. This makes this step almost instantaneous, even for large cohorts.


## Set up your environment
To be consistent in your analysis, create output folders for each of the analysis steps as follows. You can choose the change the name of these folders but you will have to change our code accordingly.
<div class="code-example" markdown="1">
```bash
dx mkdir -p Phasing/PhasingWGS/step0_qc/chunks/support
```
</div>


## Phasing the WGS data

### Quality control
We perform quality control of the variant sites and filtered out SNPs and indels for (i) Hardy-Weinberg p-value < 10-30, (ii) more than 10% of the individuals having no data (GQ score=0; missing data), (iii) heterozygous excess less than 0.5 or greater than 1.5, and (iv) alternative alleles with AAscore < 0.5. Additionally, we keep only variant sites with the tag "FILTER=PASS", as suggested by the data providers (*Halldorsson et al., Nature 2022*).

The following code is an example of our QC for the chromosomes 20 in 1288 chunks. This number has to be adapted per chromosome.




<div class="code-example" markdown="1">
```bash
# step1. Perform QC for each chunk of the chromosome
for CHR in 20; do
	JOBID0=$(dx run app-swiss-army-knife -icmd="tabix /mnt/project/Bulk/Whole\ genome\ sequences/Whole\ genome\ GraphTyper\ joint\ call\ pVCF/QC/qc_metrics_graphtyper_v2.7.1_qc.tab.gz chr${CHR} | awk 'BEGIN{print \"##fileformat=VCFv4.2\n##contig=<ID=chr${CHR}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\"}(10^$14<1 && 10^$14>0){AF=10^\$14;MAF=(AF>0.5?1-AF:AF);ExcHet=\$17/(2*(\$16+\$17+\$18)*(AF)*(1-AF));}{if (\$7>0.5 && \$15 <15012 && ExcHet >=0.5 && ExcHet <= 1.5){ split(\$3, array, \":\"); print \$1\"\t\"\$2\"\t\"\$4\"\t\"array[3]\"\t\"array[4]\"\t.\t.\t.\t.\"}}' | bcftools view -Ob -o pass_qc_chr${CHR}.bcf && bcftools index --threads 2 -f pass_qc_chr${CHR}.bcf" --tag filter_tab --tag tabix --tag swiss_army_knife --instance-type mem1_ssd1_v2_x2 --folder="Phasing/PhasingWGS/step0_qc/chunks/support" --name qc_step03a --priority normal -y | tail -n1 | cut -d" " -f3)

	#the number of chunks is manual but can easily be stored in an array using dx find data
	for CHUNK in {0..1288}; do
		START=$(((CHUNK*50000)+1))
		END=$(((CHUNK+1)*50000))
		JOBID1=$(dx run app-swiss-army-knife -iin="/Bulk/Whole\ genome\ sequences/Whole\ genome\ GraphTyper\ joint\ call\ pVCF/ukb23352_c${CHR}_b${CHUNK}_v1.vcf.gz" -iin="/Bulk/Whole\ genome\ sequences/Whole\ genome\ GraphTyper\ joint\ call\ pVCF/ukb23352_c${CHR}_b${CHUNK}_v1.vcf.gz.tbi" -iin="/data/ukb_wgs/unphased/qc/support/pass_qc_chr${CHR}.bcf" -iin="/data/ukb_wgs/unphased/qc/support/pass_qc_chr${CHR}.bcf.csi" -icmd="bcftools annotate -x ^INFO/AC,^INFO/AN,^FORMAT/GT -Ou ukb23352_c${CHR}_b${CHUNK}_v1.vcf.gz | bcftools view -f PASS -Ou | bcftools norm -m -any -Ou | bcftools view -i 'ALT!=\"*\"' -Ob -o split_ukb23352_c${CHR}_v1.bcf && bcftools index -f split_ukb23352_c${CHR}_v1.bcf && bcftools isec -c none -n=2 -w1 -r chr${CHR}:${START}-${END} split_ukb23352_c${CHR}_v1.bcf pass_qc_chr${CHR}.bcf -Ob -o pass_qc_ukb23352_c${CHR}_${START}_${END}_v1.bcf && bcftools index -f pass_qc_ukb23352_c${CHR}_${START}_${END}_v1.bcf && bcftools +fill-tags -r chr${CHR}:${START}-${END} -Ou pass_qc_ukb23352_c${CHR}_${START}_${END}_v1.bcf -- -t HWE -S /mnt/project/data/ukb_wgs/unphased/qc/support/UKB_samples_with_WGS_british_single_gbr.txt | bcftools view -G -e \"INFO/HWE_GBR < 1e-30\" -Ob -o pass_hwe_chr${CHR}_${START}_${END}.bcf && bcftools index -f pass_hwe_chr${CHR}_${START}_${END}.bcf && bcftools isec -c none -n=2 -w1 -r chr${CHR}:${START}-${END} pass_qc_ukb23352_c${CHR}_${START}_${END}_v1.bcf pass_hwe_chr${CHR}_${START}_${END}.bcf -Ob -o ukb23352_c${CHR}_${START}_${END}_v1.bcf && bcftools index -f ukb23352_c${CHR}_${START}_${END}_v1.bcf && rm -f pass* split_*" --tag qc --tag filter_tab --tag filter_hwe --instance-type mem1_ssd1_v2_x2 --folder="Phasing/PhasingWGS/step0_qc/chunks" --name qc_chr${CHR}_chunk${CHUNK} --priority normal -y | tail -n1 | cut -d" " -f3)
	done
done


# step2. Concat all chunks in a single output file
for CHR in 20; do
	dx find data --folder "Phasing/PhasingWGS/step0_qc/chunks/" --name "*c${CHR}*.bcf" --delim | sort -k 4 -t$'\t' -V | awk '{print "/mnt/project"$6}' > concat_chr${CHR}.txt
	dx upload concat_chr${CHR}.txt --path="Phasing/PhasingWGS/step0_qc/support/";
	JOBID0=$(dx run app-swiss-army-knife -icmd="bcftools concat -f /mnt/project/Phasing/PhasingWGS/step0_qc/support/concat_chr${CHR}.txt -n -o ukb23352_c${CHR}_qc_v1.bcf --threads 2 && bcftools index ukb23352_c${CHR}_qc_v1.bcf --threads 2 -f" --instance-type mem1_ssd1_v2_x2 --folder="Phasing/PhasingWGS/step0_qc/" --name concat_qc --priority normal -y | tail -n1 | cut -d" " -f3)
	JOBID1=$(dx run app-swiss-army-knife -icmd="bcftools view -G /mnt/project/Phasing/PhasingWGS/step0_qc/ukb23352_c${CHR}_qc_v1.bcf --threads 4 -Ob -o ukb23352_c${CHR}_qc_v1_sites.bcf && bcftools index ukb23352_c${CHR}_qc_v1_sites.bcf -f --threads 4" --tag sites --instance-type mem1_ssd1_v2_x4 --folder="Phasing/PhasingWGS/step0_qc/" --name get_sites --priority normal --depends-on ${JOBID0} -y | tail -n1 | cut -d" " -f3)
done
```
</div>


### Phasing
SHAPEIT5 phases common variants using the SHAPEIT5_phase_common tool. As an input, **SHAPEIT5_phase_common** requires an unphased file (with AC and AN tags), and automatically sub-sets the file to the desired MAF (e.g. 0.001). There are different strategies to phase common variants. The first, is to phase the whole chromosome in a single job. This is feasible for SNP array data and WES data, but it is not optimal for WGS data. Therefore we recommend to chunk the chromosome into large chunks (e.g. 20 cM) if using large WGS data. In the following we see how to phase WGS data in chunks.

#### Phasing common variants in chunks
When using WGS data on large sample size, it is good practice to run **SHAPEIT5_phase_common** in different large regions of the chromosomes (e.g. 20cM). In the following, we perform phasing in chunks with overlapping regions that are large enough to have a good amount of heterozygous sites for the ligation step (i.e, assembling all chunks together).


<div class="code-example" markdown="1">
```bash


```
</div>



One advantage of this approach is that these two chunks can run in parallel and they require less computational resources than the single job. However, as we wish to perform phasing at rare variants using a haplotype scaffold, we need to ligate these chunks to create a single file for the entire region.


#### Ligate chunks

#### Phasing rare variants




















































