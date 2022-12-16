---
layout: default
title: UK Biobank WGS data
nav_order: 3
parent: Tutorials
---
# UK Biobank WGS data
{: .no_toc }
<br>

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---


## Rationale
SHAPEIT5 is a two-step approach that treats each chromosome independently and works as follows:

1. Phase common variants (MAF >= 0.1%) of a chromosome using SHAPEIT5_phase_common. This can be done as a single job for SNP array or WES data, but it might be necessary to split the chromosome into large chunks (e.g. 20 cM) in the case of WGS data.

2. Ligate the phased common variants (MAF >= 0.1%) of a chromosome using SHAPEIT5_ligate, only if chunking was performed in step 1. The ligation step is computationally light and uses variants in the intersection of the chunks to provide chromosome-wide haplotypes. The result of this step (or the previous step if no chunking was used), is used as a haplotype scaffold for the next step.

3. Phase rare variants (MAF < 0.1%) of a chromosome using SHAPEIT5_rare. To do this, we use the haplotype scaffold generated in step 2 (or 1) and we proceed in relatively small chunks (e.g. 5Mb) to run relatively fast job in parallel. At the end of this step, we have several fully phased chunks across the chromosome.

4. Concatenate the phased chunks generated in step 3 using bcftools concat -n. As in the previous step haplotypes have been phased onto a haplotype scaffold, there is no need to ligate the chunks, and the files can be concatenated without decompression and recompression. This makes this step almost instantaneous, even for large cohorts.

<br>
## Phasing the WGS data
<br>
### Set up your environment
To structure the outputs of the analysis and to simplify the understanding of each step of this tutorial, we first create output folders as follows. You can choose the change the name of these folders but you will have to change the code accordingly.
<div class="code-example" markdown="1">
```bash
dx mkdir -p Phasing/PhasingWGS/step0_qc/chunks/support
dx mkdir -p Phasing/PhasingWGS/step1_phase_common/chunks
dx mkdir -p Phasing/PhasingWGS/step2_phase_rare/chunks
```
</div>
<br>

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

<br>
### Phasing
SHAPEIT5 phases common variants using the SHAPEIT5_phase_common tool. As an input, **SHAPEIT5_phase_common** requires an unphased file (with AC and AN tags), and automatically sub-sets the file to the desired MAF (e.g. 0.001). There are different strategies to phase common variants. The first, is to phase the whole chromosome in a single job. This is feasible for SNP array data and WES data, but it is not optimal for WGS data. Therefore we recommend to chunk the chromosome into large chunks (e.g. 20 cM) if using large WGS data. In the following we see how to phase WGS data in chunks.
<br>
#### Phasing common variants in chunks
When using WGS data on large sample size, it is good practice to run **SHAPEIT5_phase_common** in different large regions of the chromosomes (e.g. 20cM). In the following, we perform phasing in chunks with overlapping regions that are large enough to have a good amount of heterozygous sites for the ligation step (i.e, assembling all chunks together).

The chunks can be found here (link available soon).


**IMPORTANT**: in the following code make sure to change the shapeit5 docker image name (here `shapeit5_beta.tar.gz`) to the latest version that you've downloaded [here](https://odelaneau.github.io/shapeit5/docs/installation/docker)


<div class="code-example" markdown="1">
```bash
# step1. Download map files
dx mkdir -p data/shapeit_maps/
wget https://github.com/odelaneau/shapeit5/raw/main/maps/genetic_maps.b38.tar.gz
tar -xvzf genetic_maps.b38.tar.gz
dx upload *.b38.gmap.gz --path="data/shapeit_maps/"


# step2. Download chunk file


# step3. Phasing
for CHR in 20; do
	MAP=/mnt/project/data/shapeit_maps/chr${CHR}.b38.gmap.gz
	while read LINE; do
		REG=$(echo $LINE | awk '{ print $3; }')
		BCF=/mnt/project/Phasing/PhasingWGS/step0_qc/ukb23352_c${CHR}_qc_v1.bcf
		
		OUT=ukb23352_c${CHR}_qc_v1.$REG.shapeit5.bcf
		LOG=ukb23352_c${CHR}_qc_v1.$REG.shapeit5.log
		TIM=ukb23352_c${CHR}_qc_v1.$REG.shapeit5.time
		
		dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_beta.tar.gz" --folder="/Phasing/PhasingWGS/step1_phase_common/chunks/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase_common --input $BCF --map $MAP --output $OUT --thread 32 --log $LOG --filter-maf 0.001 --region $REG && bcftools index -f $OUT --threads 32" --instance-type mem3_ssd1_v2_x32 --priority normal --name WGS_shapeit5_common_$REG -y
		
	done < ${CHUNK_FILE}
done

```
</div>

The full list of options for the **SHAPEIT5_phase_common** command can be found [here](https://odelaneau.github.io/shapeit5/docs/documentation/phase_common/).


One advantage of this approach is that these two chunks can run in parallel and they require less computational resources than the single job. However, as we wish to perform phasing at rare variants using a haplotype scaffold, we need to ligate these chunks to create a single file for the entire region.

<br>
#### Ligate chunks
Ligation of the phased common variants is performed using the SHAPEIT5_ligate tool. The program requires an ordered list of the phased common variants. We recommend using appropriate naming of the files in the previous step, so that a command such as ls -1v can directly produce the list of files in the right order.


<div class="code-example" markdown="1">
```bash
CHR=20
#generate list for SHAPEIT5_ligate
dx ls -1v /Phasing/PhasingWGS/step1_phase_common/chunks/ukb23352_c${CHR}_qc_v1.*.shapeit5.bcf > list_ligate.chr${CHR}.txt
dx upload list_ligate.chr${CHR}.txt --path="/Phasing/PhasingWGS/step1_phase_common/"

#ligate the chunks
dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_beta.tar.gz" --folder="/Phasing/PhasingWGS/step1_phase_common/" -icmd="SHAPEIT5_ligate --input /mnt/proejct/Phasing/PhasingWGS/step1_phase_common/list_ligate.chr${CHR}.txt --output ukb23352_c${CHR}_qc_v1.bcf --thread 32 --index" --instance-type mem3_ssd1_v2_x32 --priority normal --name WGS_shapeit5_chr${CHR}_ligate -y
```
</div>

At the end of this step, we have a region-wide haplotype scaffold.

<br>
#### Phasing rare variants in chunks
Having obtained the haplotype scaffold at common variants in the previous step, we can now phase rare variants using the **SHAPEIT5_phase_rare** tool. We always do this in chunks, in order to maximise the parallelization of our jobs (e.g. in 5Mb chunks). The program requires two input files, one is the haplotype scaffold generated in the previous steps, and the other input required is a file containing the whole unphased region, the same input file provided for **SHAPEIT5_phase_common**: SHAPEIT5 automatically discards common variants for this file, duplicated in the phased scaffold.




		
<div class="code-example" markdown="1">
```bash
for CHR in 20; do
	MAP=/mnt/project/data/shapeit_maps/chr${CHR}.b38.gmap.gz
	while read LINE; do

		BCF=/mnt/project/Phasing/PhasingWGS/step0_qc/ukb23352_c${CHR}_qc_v1.bcf
		SCAF=/mnt/project/Phasing/PhasingWGS/step1_phase_common/ukb23352_c20_qc_v1.bcf
		

		SCAFFOLD_REG=$(echo $LINE | awk '{ print $3; }')
		SCAFFOLD_REG_START=$(echo ${SCAFFOLD_REG} | cut -d":" -f 2 | cut -d"-" -f1)
		SCAFFOLD_REG_END=$(echo ${SCAFFOLD_REG} | cut -d":" -f 2 | cut -d"-" -f2)
		SCAFFOLD_REG_NAME=${CHR}_${SCAFFOLD_REG_START}_${SCAFFOLD_REG_END}
			
		INPUT_REG=$(echo $LINE | awk '{ print $4; }')
		INPUT_REG_START=$(echo ${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f1)
		INPUT_REG_END=$(echo ${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f2)
		INPUT_REG_NAME=${CHR}_${INPUT_REG_START}_${INPUT_REG_END}


		OUT=ukb23352_c${CHR}_qc_v1.$REG.shapeit5.bcf
		LOG=ukb23352_c${CHR}_qc_v1.$REG.shapeit5.log
		TIM=ukb23352_c${CHR}_qc_v1.$REG.shapeit5.time
		
		dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_beta.tar.gz" --folder="/Phasing/PhasingWGS/step2_phase_rare/chunks/" -icmd="/usr/bin/time -vo $TIM SHAPEIT5_phase_rare --input-plain $BCF --map $MAP --scaffold $SCAF --output $OUT --thread 32 --log $LOG --filter-maf 0.001 --scaffold-region $SCAFFOLD_REG --input-region $INPUT_REG && bcftools index -f $OUT --threads 32" --instance-type mem3_ssd1_v2_x32 --priority normal --name WGS_shapeit5_rare_$REG -y
		
	done < ${CHUNK_FILE}
done

```
</div>	

The full list of options for the **SHAPEIT5_phase_rare** command can be found [here](https://odelaneau.github.io/shapeit5/docs/documentation/phase_rare/).

<br>

#### Concatenate chunks
Since the phasing of rare variants has been performed using a ligated haplotype scaffold, the chunks can now be concatenated together to resolve the phasing for the entire chromosome. This is done using the **bcftools concat** command.


<div class="code-example" markdown="1">
```bash
CHR=20
#generate list for bcftools concat
dx ls -1v /Phasing/PhasingWGS/step2_phase_rare/chunks/ukb23352_c${CHR}_qc_v1.*.shapeit5.bcf > list_ligate.chr${CHR}.txt
dx upload list_ligate.chr${CHR}.txt --path="/Phasing/PhasingWGS/step1_phase_rares/"

#concatenate chunks
OUT=ukb23352_c${CHR}_qc_v1.phased.bcf
dx run app-swiss-army-knife --folder="/Phasing/PhasingWGS/step2_phase_rare/" -icmd="bcftools concat -n -f /mnt/project/Phasing/PhasingWGS/step1_phase_rares/list_ligate.chr${CHR}.txt -Ob -o $OUT --threads 32 && bcftools index -f $OUT --threads 32" --instance-type mem1_ssd1_v2_x32 --priority normal --name WGS_shapeit5_rare_chr${CHR}_concat -y	

```
</div>	



<br>


## Validation of your phasing
You can validate the quality of the haplotypes using the **SHAPEIT5_switch** tool. For this you will need parent-offspring duos or trios stored in a three-columns file (here called `family.ped`) following this format:

<div class="code-example" markdown="1">
- family.ped:      `offspring_id`    `parent1_id`    `parent2_id`

</div>


To validate the phasing using family data, the phasing must be performed by excluding parental genomes, so that offsprings are phased regardless of their parental genomes. This can be done using the **bcftools view** command, with as input the output data of the quality control step (located here `/mnt/project/Phasing/PhasingWGS/step0_qc/ukb23352_c${CHR}_qc_v1.bcf`).

<div class="code-example" markdown="1">
```bash
dx mkdir -p Phasing/PhasingWGS/benchmark/
for CHR in {1..22}; do
IN=/mnt/project/Phasing/PhasingWGS/step0_qc/ukb23352_c${CHR}_qc_v1.bcf
OUT=benchmark_ukb23352_c${CHR}_qc_v1.bcf
dx run app-swiss-army-knife --folder "/Phasing/PhasingWGS/benchmark/" -icmd="bcftools view --threads 16 -S ^parents.txt -Ob -o ${OUT} ${IN} && bcftools index ${OUT} --threads 16" --instance-type mem1_ssd1_v2_x16 --priority normal --name phasing_chr${CHR} -y
done
```
</div>

After excluding parental genomes with the above command, proceed with the normal phasing procedure as detailed above.

Let's consider that you performed the above steps of phasing using the input data exluding parental genomes, which produced a phased output file that you named `benchmark_ukb23352_c${CHR}_qc_v1.phased.bcf`). You can validate you phasing using the following command:


<div class="code-example" markdown="1">
```bash
for CHR in 20; do
dx run app-swiss-army-knife --folder "/Phasing/PhasingWGS/benchmark/" -iimage_file="/docker/shapeit5_beta.tar.gz" -icmd="SHAPEIT5_switch --validation /mnt/project/Phasing/PhasingWGS/step0_qc/ukb23352_c${CHR}_qc_v1.bcf --estimation benchmark_ukb23352_c${CHR}_qc_v1.phased.bcf --region ${CHR} --output benchmark_chr${CHR}"  --instance-type mem2_ssd1_v2_x16 --priority normal --name benchmark_chr${CHR} -y
done
```
</div>



The full list of options for the **SHAPEIT5_switch** command can be found [here](https://odelaneau.github.io/shapeit5/docs/documentation/switch/).



















