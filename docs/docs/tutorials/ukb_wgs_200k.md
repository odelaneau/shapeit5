---
layout: default
title: UK Biobank WGS data - 200k release
nav_order: 3
parent: Tutorials
---
# UK Biobank WGS data
{: .no_toc }


## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

**For any question on this pipeline, please contact [Robin J. Hofmeister](https://github.com/RJHFMSTR) and [Diogo M. Ribeiro](https://diogomribeiro.github.io/).** 

---

## Rationale
This tutorial provides a brief explanation of the scripts that were used to estimate the haplotype of the newest release of the UK Biobank WGS data, interim release of 200,031 individuals. A phasing report documenting the quality control steps, a more detailed rationale and the validation of the phasing is available [here](https://docs.google.com/document/d/1EJmh-JcR8HBvu3rjBtREw_50kDIW-sOW2EcC6zweuTc/edit?usp=sharing).

<br>
## Phasing the WGS data - interim release of 200,031 individuals
<br>

SHAPEIT5 is a two-step approach that treats each chromosome independently and works as follows:

1. Phase common variants (MAF >= 0.1%) of a chromosome using SHAPEIT5_phase_common. This can be done as a single job for SNP array or WES data, but it might be necessary to split the chromosome into large chunks (e.g. 20 cM) in the case of WGS data. Chunk coordinates for b38 are given in the resources folder.

2. Ligate the phased common variants (MAF >= 0.1%) of a chromosome using SHAPEIT5_ligate, only if chunking was performed in step 1. The ligation step is computationally light and uses variants in the intersection of the chunks to provide chromosome-wide haplotypes. The result of this step (or the previous step if no chunking was used), is used as a haplotype scaffold for the next step.

3. Phase rare variants (MAF < 0.1%) of a chromosome using SHAPEIT5_rare. To do this, we use the haplotype scaffold generated in step 2 (or 1) and we proceed in relatively small chunks (e.g. 5Mb) to run relatively fast job in parallel. At the end of this step, we have several fully phased chunks across the chromosome.

4. Concatenate the phased chunks generated in step 3 using bcftools concat -n. As in the previous step haplotypes have been phased onto a haplotype scaffold, there is no need to ligate the chunks, and the files can be concatenated without decompression and recompression. This makes this step almost instantaneous, even for large cohorts.



In addition to what has been described in the previous tutorials (i.e phasing UKB WES and UKB WES - interim 150k release), this new pipeline include two novelties:

- **family phasing** : we leveraged the cryptic relatedness of the UK Biobank data to phase offspring from parental genome. This is applied chromosome-wide, but also genome-wide. It means that across all chromosomes, the first haplotypes are always inherited from the first parent listed in the pedigree file. See the [documentation](https://odelaneau.github.io/shapeit5/docs/documentation/phase_common/#usage2-phasing-related-samples) for more informations.
- **variable ploidy** : Phasing process accounts for haploidy of males on chromosome X.

<br>
It is important to note that a genome-wide chunking has been produce for this phasing. We used large chunks (~20-25MB) to phase common variants and smaller chunks (~4-7MB) to phase rare variants. This chunking has been optimized for phasing accuracy and computational costs. We based our chunking algorithm on the number of variant per chunk, the size of chunks in megabase and the size of chunks in centimorgan, to have similar chunk parameters along the genome. It can easyly be used on another dataset in genome built GRCh38. The chunking can be found in our github repository by navigating into the `chunks/` folder.

<br>


For a more detail description of the steps, see below.


<br>
<br>
## Phasing autosomes

### Quality control
We perform quality control of the variant sites and filtered out SNPs and indels for (i) Hardy-Weinberg p-value < 10-30, (ii) more than 10% of the individuals having no data (GQ score=0; missing data), (iii) alternative allele undetermined, and (iv) alternative alleles with AAscore < 0.8. Additionally, we keep only variant sites with the tag "FILTER=PASS", as suggested by the data providers (*Halldorsson et al., Nature 2022*).

The following code was used to run the QC on all UK Biobank WGS chunks.

The data is initially provided in more than 60k chunks. We first grouped these chunks in ~100 batches per chromosomes. This allows us to run an entire chromosome QC in a hundred jobs in parallel on the UK Biobank RAP. Whitin each of these 100 jobs, we used the xargs command to parallelize 8 jobs at the same time. This code is provided below. It relies on (i) a QC script (can be found below or [here](https://github.com/odelaneau/shapeit5/tree/main/tasks/phasingUKB_200k_release/autosomes)), (ii) QC chunks (can be found [here](https://github.com/odelaneau/shapeit5/tree/main/tasks/phasingUKB_200k_release/autosomes)).


<div class="code-example" markdown="1">
```bash
#!bin/bash


CHR_START=$1
CHR_END=$2



for CHR in $(seq $CHR_START 1 $CHR_END); do

	ODIR0=QC/chr${CHR}/support
	dx mkdir -p ${ODIR0}


	SDIR=/QC/support
	# 1. upload script
	SCR=script_stats_qc.sh # This script is provided below. Can also be found on the SHAPEIT5 github /tasks/phasingUKB_200k_release/autosomes
	dx rm ${SDIR}/${SCR}
	dx upload ${SCR} --path="${SDIR}/"
	
	# 2. upload sample file
	SAMP=ukb_200k.samples_with_pop_tag.unrelated_only.txt # this is a two column file with col1= sample_id, col2=population. It contains only unrelated individuals. For col2, we used CAU for caucasians et UN for the remaining samples. It allows us to use caucasian metric only for QC (e.g HWE).

	

	N=$(($(ls chunks_qc/Splitted/chr${CHR}/chr${CHR}_chunks* | wc -l)-1)) # the data before QC is stored in more thank 60k small chunks genome-wide. We groupped these in 100 groups of chunks for each chromosome (since we can run ~100 jobs in parrallel on the RAP). We run on QC job per group of chunks. Within each of these jobs, 8 QC jobs run in parallel (using xargs command). This allows us to considerably speed up the genome-wide running time.


	for SPL in $(seq 0 1 ${N}); do
		echo ${CHR}; echo ${SPL}


		printf -v ID "%03d" $(echo $SPL)
	
		# 3. upload chunk file
		PARAM=chunks_qc/Splitted/chr${CHR}/chr${CHR}_chunks${ID}
#		dx rm "${ODIR0}/$(basename ${PARAM})"
		dx upload ${PARAM} --path="${ODIR0}/"
		PARAM=$(basename ${PARAM})
	
		# 4. run script in parallele
		ODIR1=QC/chr${CHR}/chunks
		dx mkdir -p ${ODIR1}
	
		# 4.1 concat the chunks.
		IN1=ukb24304_chr${CHR}_b*_v1.qceed.bcf
		OUT1=ukb24304_chr${CHR}_chunk_${SPL}.qceed.bcf
		IN2=ukb24304_chr${CHR}_b*_v1.NOqceed.stats.bcf
		OUT2=ukb24304_chr${CHR}_chunk_${SPL}.NOqceed.stats.bcf
		threads=16
		
		dx run app-swiss-army-knife -iin="${ODIR0}/${PARAM}" -iin="${SDIR}/${SCR}" -iin="${SDIR}/${SAMP}" -icmd="cat ${PARAM} | xargs -P 8 -n2 bash ${SCR} && ls -1v ${IN1} > concat_list_chr${CHR}.txt && bcftools concat --threads ${threads} --naive-force -f concat_list_chr${CHR}.txt -o ${OUT1} && bcftools index ${OUT1} --threads ${threads} && rm concat_list_chr${CHR}.txt && rm ${IN1}* && ls -1v ${IN2} > concat_list2_chr${CHR}.txt && bcftools concat --threads ${threads} --naive-force -f concat_list2_chr${CHR}.txt -o ${OUT2} && bcftools index ${OUT2} --threads ${threads} && rm concat_list2_chr${CHR}.txt && rm ${IN2}*" --instance-type mem2_ssd1_v2_x16 --folder="${ODIR1}" --name QC_chr${CHR}_spl_${SPL} --priority normal -y --brief --ignore-reuse

	
	done
done

```
</div>
<br>

The following code corresponds to the script named `script_stats_qc.sh` above. Can also be found [here](https://github.com/odelaneau/shapeit5/tree/main/tasks/phasingUKB_200k_release/autosomes).

<div class="code-example" markdown="1">
```bash
#!bin/bash

CHR=$1
CNK=$2

echo $CHR; echo $CNK

IN=/mnt/project/Bulk/Whole\ genome\ sequences/Population\ level\ WGS\ variants\,\ pVCF\ format\ -\ interim\ 200k\ release/ukb24304_c${CHR}_b${CNK}_v1.vcf.gz
SAMP=ukb_200k.samples_with_pop_tag.unrelated_only.txt
OUT=ukb24304_chr${CHR}_b${CNK}_v1.qceed.bcf
OUT2=ukb24304_chr${CHR}_b${CNK}_v1.NOqceed.stats.bcf
TMP=ukb24304_chr${CHR}_b${CNK}_v1.tmp.bcf
threads=2


bcftools norm --threads ${threads} -m -any "${IN}" -Ou | bcftools annotate -x "FORMAT" -Ou | bcftools +fill-tags --threads ${threads} -Ob -o ${TMP} -- -t HWE,AF,ExcHet -S ${SAMP} && bcftools view ${TMP} --threads ${threads} -f PASS -e 'ALT="*" | F_MISSING>0.1 | INFO/AAScore<0.8 | INFO/AF=0 | INFO/AF=1 | INFO/HWE_CAU<1e-30' -Ob -o ${OUT} && bcftools index ${OUT} --threads ${threads} && bcftools view -G ${TMP} --threads ${threads} -Ob -o ${OUT2} && bcftools index ${OUT2} --threads ${threads} && rm ${TMP}


```
</div>




<br>
### Phasing
SHAPEIT5 phases common variants using the SHAPEIT5_phase_common tool. As an input, **SHAPEIT5_phase_common** requires an unphased file (with AC and AN tags), and automatically sub-sets the file to the desired MAF (e.g. 0.001). There are different strategies to phase common variants. The first, is to phase the whole chromosome in a single job. This is feasible for SNP array data and WES data, but it is not optimal for WGS data. Therefore we recommend to chunk the chromosome into large chunks (e.g. 20 cM) if using large WGS data. In the following we see how to phase WGS data in chunks.
<br>
#### Phasing common variants in chunks
When using WGS data on large sample size, it is good practice to run **SHAPEIT5_phase_common** in different large regions of the chromosomes (e.g. 20cM). In the following, we perform phasing in chunks with overlapping regions that are large enough to have a good amount of heterozygous sites for the ligation step (i.e, assembling all chunks together).

The chunks can be found [here](https://github.com/odelaneau/shapeit5/tree/main/tasks/phasingUKB_200k_release/autosomes).

We used the option `--pedigree` to leverage family information (i.e parent-offspring duos and trios) in the phasing process. See the [documentation](https://odelaneau.github.io/shapeit5/docs/documentation/phase_common/#usage2-phasing-related-samples) for more informations.


**IMPORTANT**: in the following code make sure to change the shapeit5 docker image name (here `shapeit5.ukb200k.tar.gz`) to the latest version that you've downloaded [here](https://odelaneau.github.io/shapeit5/docs/installation/docker)


<div class="code-example" markdown="1">
```bash
#!/bin/bash

CHR_START=$1
CHR_END=$2


PREFIX=ukb200k

BIN=docker/shapeit5.ukb200k.tar.gz # if not done yet, upload the SHAPEIT5 docker image onto the RAP.


# step0. Create output directory.
ODIR=Phasing/step1_phase_common/chunks
dx mkdir -p ${ODIR}

# step1. Upload map files
S5_DIR=../Git_repository/shapeit5 ## specify your own path for the SHAPEIT5 github, or download the map files separately.
cp ${S5_DIR}/maps/genetic_maps.b38.tar.gz ./
tar xvzf genetic_maps.b38.tar.gz
dx mkdir -p data/shapeit_maps/
dx upload *.b38.gmap.gz --path="data/shapeit_maps/"


PED=/mnt/project/data/ukb_200k.ped # the pedifree file contains three columns: col1=offspring; col2=father; col3=mother. NA can be specified for unknown parent. Parent-offspring trios and duos are determined from the UK Biobank relatedness file (provided as part of the official released - computed using the KING software). We used 0.1717<Kinship<0.3535 & IBS0<0.0012.


# step2. phasing common variants (MAF>0.001)
for CHR in $(seq $CHR_START 1 $CHR_END); do


	CHUNKS=chunks/phase_common/chunks_chr${CHR}.txt
	BCF=/mnt/project/QC/chr${CHR}/ukb24304_chr${CHR}.qceed.bcf # QCeed data, unphased, full chromosome, full sample set.
	MAP=/mnt/project/data/shapeit_maps/chr${CHR}.b38.gmap.gz

	while read LINE; do
		REG=$(echo $LINE | awk '{ print $3; }')
		CHUNK_NBR=$(echo $LINE | awk '{ print $1; }')
		OUT=${PREFIX}_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_common.bcf
		LOG=${PREFIX}_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_common.log
		TIM=${PREFIX}_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_common.time
		
		dx run app-swiss-army-knife -iimage_file="/${BIN}" --folder="${ODIR}/" -icmd="/usr/bin/time -vo $TIM phase_common_static --input $BCF --map $MAP --output $OUT --thread 72 --log $LOG --filter-maf 0.001 --region $REG --pedigree ${PED} && bcftools index -f $OUT --threads 72" --instance-type mem1_ssd1_v2_x72 --priority normal --name ${PREFIX}_shapeit5_common_chr${CHR}_${CHUNK_NBR} -y --brief --ignore-reuse
		
	done < ${CHUNKS}
done


```
</div>


<br>
#### Ligate chunks
Ligation of the phased common variants is performed using the SHAPEIT5_ligate tool. The program requires an ordered list of the phased common variants. We recommend using appropriate naming of the files in the previous step, so that a command such as ls -1v can directly produce the list of files in the right order.


<div class="code-example" markdown="1">
```bash
#!bin/bash


PREFIX=ukb200k

mkdir -p TMP

BIN=docker/shapeit5.ukb200k.tar.gz
IDIR=Phasing/step1_phase_common/chunks
ODIR=Phasing/step2_ligate
dx mkdir -p ${ODIR}
threads=8

PED=ukb_200k.ped # sample file as in step2_phase_common.sh

for CHR in 7; do

	SCR=TMP/script_ligate_chr${CHR}.sh
	IN=/mnt/project/${IDIR}/UKB_chr${CHR}.chunk_*.shapeit5_common.bcf
	OUT=${PREFIX}_chr${CHR}.shapeit5_common_ligate.bcf


	printf '#!bin/bash\ndx run app-swiss-army-knife --folder='${ODIR}/' -iimage_file='/${BIN}' ' > ${SCR}
	
	N=$(($(cat chunks/phase_common/chunks_chr${CHR}.txt | wc -l)-1))

	LINE=''

	for SPL in $(seq 0 1 ${N}); do

		LINE="$LINE"-iin="${IDIR}"/"${PREFIX}"_chr"${CHR}".chunk_"${SPL}".shapeit5_common.bcf" "
		LINE="$LINE"-iin="${IDIR}"/"${PREFIX}"_chr"${CHR}".chunk_"${SPL}".shapeit5_common.bcf.csi" "
		
	done

	LINE2='-iin=data/'${PED}' -icmd="ls -1v '${PREFIX}'_chr'${CHR}'.chunk_*.shapeit5_common.bcf > list_ligate.chr'${CHR}'.txt && ligate_static --pedigree '${PED}' --input list_ligate.chr'${CHR}'.txt --output '${OUT}' --thread '${threads}' --index" --instance-type mem1_ssd1_v2_x8 --priority low --name ligate_chr'${CHR}' -y --brief --ignore-reuse'

	L="$LINE""$LINE2"
	echo ${L} >> ${SCR}
	bash ${SCR}

done

```
</div>

At the end of this step, we have a region-wide haplotype scaffold.

<br>
#### Phasing rare variants in chunks
Having obtained the haplotype scaffold at common variants in the previous step, we can now phase rare variants using the **SHAPEIT5_phase_rare** tool. We always do this in chunks, in order to maximise the parallelization of our jobs (e.g. in 5Mb chunks). The program requires two input files, one is the haplotype scaffold generated in the previous steps, and the other input required is a file containing the whole unphased region, the same input file provided for **SHAPEIT5_phase_common**: SHAPEIT5 automatically discards common variants for this file, duplicated in the phased scaffold.

To minimize running times and costs, we grouped the phasing chunks into ~70 batches of 8 chunks. We ran each batch as a single job on the UK Biobank RAP. Whithin each job, we used the xargs command to parallelize 4 jobs running at the same time on the same machine.


		
<div class="code-example" markdown="1">
```bash
#!/bin/bash

# Split the total number of phase_rare chunks to run several jobs in parallel on the same machine.
n=8
mkdir -p chunks_rare/Splitted
rm chunks_rare/Splitted/*

cp chunks/phase_rare/chunks_chr*.txt chunks_rare/
cat chunks_rare/chunks_chr*.txt | cut -f -4 > chunks_rare/chunks.txt
split -d -a 3 -l ${n} chunks_rare/chunks.txt chunks_rare/Splitted/chunk_
N=$(($(ls chunks_rare/Splitted/chunk_* | wc -l)-1))


# phasing rare variants: 4 chunks_rare in parallel, 8 threads per job, 32 cpu total

PREFIX=ukb200k
BIN=docker/shapeit5.ukb200k.tar.gz
threads=8
IDIR=Phasing/step2_ligate
SDIR=Phasing/step3_phase_rare/support
ODIR=Phasing/step3_phase_rare/chunks
dx mkdir -p ${ODIR}
dx mkdir -p ${SDIR}
PED=ukb_200k.ped


SCR=script_phase_rare.sh
dx upload ${SCR} --path="${SDIR}/"



for SPL in $(seq 1 1 ${N}); do

	printf -v ID "%03d" $(echo $SPL)

	PARAM=chunks_rare/Splitted/chunk_${ID}
	dx rm "${SDIR}/$(basename ${PARAM})"
	dx upload ${PARAM} --path="${SDIR}/"
	PARAM=$(basename ${PARAM})
		
	dx run app-swiss-army-knife -iimage_file="/${BIN}" -iin="${SDIR}/${PARAM}" -iin="${SDIR}/${SCR}" -iin="data/${PED}" -icmd="cat ${PARAM} | xargs -P 4 -n4 bash ${SCR}" --instance-type mem3_ssd1_v2_x32 --folder="${ODIR}" --name shapeit5_rare_chunk_${SPL} --priority normal -y --brief --ignore-reuse


done

```
</div>	

The full list of options for the **SHAPEIT5_phase_rare** command can be found [here](https://odelaneau.github.io/shapeit5/docs/documentation/phase_rare/).

<br>

#### Concatenate chunks
Since the phasing of rare variants has been performed using a ligated haplotype scaffold, the chunks can now be concatenated together to resolve the phasing for the entire chromosome. This is done using the **bcftools concat** command.


<div class="code-example" markdown="1">
```bash
#!/bin/bash

PREFIX=ukb200k

threads=16

IDIR=Phasing/step3_phase_rare/chunks
ODIR=Phasing/step3_phase_rare

for CHR in {1..22}; do
	chunks=/mnt/project/${IDIR}/${PREFIX}_chr${CHR}.chunk_*.shapeit5_rare.bcf
	OUT=${PREFIX}_chr${CHR}.full.shapeit5_rare.bcf
	
	dx run app-swiss-army-knife --folder="${ODIR}/" -icmd="ls -1v ${chunks} > concat_list_chr${CHR}.txt && bcftools concat -n -f concat_list_chr${CHR}.txt -o ${OUT} && bcftools index ${OUT} --threads ${threads} && rm concat_list_chr${CHR}.txt" --instance-type mem1_ssd1_v2_x16 --priority normal --name WGS_concat_chr${CHR} -y
	
done


```
</div>	



<br>




## Phasing chromosome X

### Quality control
Similarly to the autosomes, we perform quality control of the variant sites. On chromosome X, we used only female for the population-based metrics, such as HWE.


<div class="code-example" markdown="1">
```bash
#!bin/bash

for CHR in X; do

	ODIR0=QC/chr${CHR}/support
	dx mkdir -p ${ODIR0}


	SDIR=QC/support
	# 1. upload script
	SCR=script_stats_qc.chrX.sh
	dx upload ${SCR} --path="${SDIR}/"
	
	# 2. upload sample file
	SAMP=ukb_200k.samples_HWE_females.txt # similar as for the autosomes, we computed population metrics on caucasians only. Specifically, on chromosome X, we used only females. This file contains two columns: col1=sample id; col2=population tag. For col2, we used CAU_F for caucasian females, CAU_M for caucasian males, UN_F for the non-caucasian females and UN_M for the non-caucasian males.

	N=$(($(ls ../autosomes/chunks_qc/Splitted/chr${CHR}/chr${CHR}_chunks* | wc -l)-1))



	for SPL in $(seq 0 1 3); do
		echo ${CHR}; echo ${SPL}


		printf -v ID "%03d" $(echo $SPL)
	
		# 3. upload chunk file
		PARAM=../chunks_qc/Splitted/chr${CHR}/chr${CHR}_chunks${ID}
		dx rm "${ODIR0}/$(basename ${PARAM})"
		dx upload ${PARAM} --path="${ODIR0}/"
		PARAM=$(basename ${PARAM})
	
		# 4. run script in parallele
		ODIR1=QC/chr${CHR}/chunks
		dx mkdir -p ${ODIR1}
	
		# 4.1 concat the chunks.
		IN1=ukb24304_chr${CHR}_b*_v1.qceed.bcf
		OUT1=ukb24304_chr${CHR}_chunk_${SPL}.qceed.bcf
		IN2=ukb24304_chr${CHR}_b*_v1.NOqceed.stats.bcf
		OUT2=ukb24304_chr${CHR}_chunk_${SPL}.NOqceed.stats.bcf
		threads=16
		
		dx run app-swiss-army-knife -iin="${ODIR0}/${PARAM}" -iin="${SDIR}/${SCR}" -iin="${SDIR}/${SAMP}" -icmd="cat ${PARAM} | xargs -P 8 -n2 bash ${SCR} && ls -1v ${IN1} > concat_list_chr${CHR}.txt && bcftools concat --threads ${threads} --naive-force -f concat_list_chr${CHR}.txt -o ${OUT1} && bcftools index ${OUT1} --threads ${threads} && rm concat_list_chr${CHR}.txt && rm ${IN1}* && ls -1v ${IN2} > concat_list2_chr${CHR}.txt && bcftools concat --threads ${threads} --naive-force -f concat_list2_chr${CHR}.txt -o ${OUT2} && bcftools index ${OUT2} --threads ${threads} && rm concat_list2_chr${CHR}.txt && rm ${IN2}*" --instance-type mem2_ssd1_v2_x16 --folder="${ODIR1}" --name QC_chr${CHR}_spl_${SPL} --priority normal -y --brief --ignore-reuse

	
	done
done

```
</div>
<br>

The following code corresponds to the script named `script_stats_qc.sh` above. Can also be found [here](https://github.com/odelaneau/shapeit5/tree/main/tasks/phasingUKB_200k_release/chrX).

<div class="code-example" markdown="1">
```bash
#!bin/bash

CHR=$1
CNK=$2

echo $CHR; echo $CNK

IN=/mnt/project/Bulk/Whole\ genome\ sequences/Population\ level\ WGS\ variants\,\ pVCF\ format\ -\ interim\ 200k\ release/ukb24304_c${CHR}_b${CNK}_v1.vcf.gz
SAMP=ukb_200k.samples_HWE_females.txt
OUT=ukb24304_chr${CHR}_b${CNK}_v1.qceed.bcf
OUT2=ukb24304_chr${CHR}_b${CNK}_v1.NOqceed.stats.bcf
TMP=ukb24304_chr${CHR}_b${CNK}_v1.tmp.bcf
threads=2


bcftools norm --threads ${threads} -m -any "${IN}" -Ou | bcftools annotate -x "FORMAT" -Ou | bcftools +fill-tags --threads ${threads} -Ob -o ${TMP} -- -t HWE,AF,ExcHet -S ${SAMP} && bcftools view ${TMP} --threads ${threads} -f PASS -e 'ALT="*" | INFO/AAScore<0.8 | INFO/AF=0 | INFO/AF=1 | INFO/HWE_CAU_F<1e-30' -Ob -o ${OUT} && bcftools index ${OUT} --threads ${threads} && bcftools view -G ${TMP} --threads ${threads} -Ob -o ${OUT2} && bcftools index ${OUT2} --threads ${threads} && rm ${TMP}


```
</div>

Additionally, since we modelled variably ploidy on chromosome X, we removed PAR regions.

<div class="code-example" markdown="1">
```bash
#!bin/bash

threads=32
ODIR=QC/chrX/PAR_filter/
dx mkdir -p ${ODIR}

# filter QCeed data.
IN=/mnt/project/QC/chrX/ukb24304_chrX.qceed.bcf
REG0=chrX:0-10000
OUT0=ukb24304_chrX_0_10000.bcf
dx run app-swiss-army-knife -icmd="bcftools view -i 'F_MISSING<0.1' --threads ${threads} -r ${REG0} -Ob -o ${OUT0} ${IN} && bcftools index ${OUT0} --threads ${threads}" --instance-type mem2_ssd1_v2_x32 --folder="${ODIR}" --name PAR0_chr${CHR} --priority normal -y --brief

REG1=chrX:2781480-155701384
OUT1=ukb24304_chrX_2781480_155701384.bcf
dx run app-swiss-army-knife -icmd="bcftools view -i 'F_MISSING<0.1' --threads ${threads} -r ${REG1} -Ob -o ${OUT1} ${IN} && bcftools index ${OUT1} --threads ${threads}" --instance-type mem2_ssd1_v2_x32 --folder="${ODIR}" --name PAR1_chr${CHR} --priority normal -y --brief

REG2=chrX:156030896-1000000000
OUT2=ukb24304_chrX_156030896_1000000000.bcf
dx run app-swiss-army-knife -icmd="bcftools view -i 'F_MISSING<0.1' --threads ${threads} -r ${REG2} -Ob -o ${OUT2} ${IN} && bcftools index ${OUT2} --threads ${threads}" --instance-type mem2_ssd1_v2_x32 --folder="${ODIR}" --name PAR2_chr${CHR} --priority normal -y --brief

# merge:
OUT3=ukb24304_chrX_without_PAR.bcf
ODIR3=QC/chrX/
dx run app-swiss-army-knife -icmd="bcftools concat -n -Ob -o ${OUT3} /mnt/project/${ODIR}/${OUT0} /mnt/project/${ODIR}/${OUT1} /mnt/project/${ODIR}/${OUT2}&& bcftools index ${OUT3} --threads ${threads}" --instance-type mem2_ssd1_v2_x32 --folder="${ODIR3}" --name PAR2_chr${CHR} --priority normal -y --brief

OUT4=ukb24304_chrX_without_PAR.sites.bcf
dx run app-swiss-army-knife -icmd="bcftools view -G --threads ${threads} -Ob -o ${OUT4} /mnt/project/${ODIR3}/${OUT3} && bcftools index ${OUT4} --threads ${threads}" --instance-type mem2_ssd1_v2_x32 --folder="${ODIR3}" --name PAR2_chr${CHR} --priority normal -y --brief
```
</div>


<br>

### Phasing
For the phasing, we proceed similarly as for autosome. The main difference is the option `--haploid` which allows to specify the list of haploid males. More details about our procedure to determine this list of individuals can be found in the [phasing report](https://docs.google.com/document/d/1EJmh-JcR8HBvu3rjBtREw_50kDIW-sOW2EcC6zweuTc/edit?usp=sharing).
<br>
#### Phasing common variants in chunks


<div class="code-example" markdown="1">
```bash
#!/bin/bash

PREFIX=ukb200k
BIN=docker/shapeit5.ukb200k.tar.gz

# step0. Create output directory.
ODIR=Phasing/step1_phase_common/chunks
dx mkdir -p ${ODIR}

# step1. Upload map files
#cp ../Git_repository/shapeit5/maps/genetic_maps.b38.tar.gz ./
#tar xvzf genetic_maps.b38.tar.gz
#dx mkdir -p data/shapeit_maps/
#dx upload *.b38.gmap.gz --path="data/shapeit_maps/"

PED=/mnt/project/data/ukb_200k.ped

# step2. phasing common variants (MAF>0.001)
for CHR in X; do

	CHUNKS=../autosomes/chunks/phase_common/chunks_chr${CHR}.txt
	BCF=/mnt/project/QC/chr${CHR}/ukb24304_chrX_without_PAR.bcf # QCeed data, unphased, full chromosome, full sample set.
	MAP=/mnt/project/data/shapeit_maps/chr${CHR}.b38.gmap.gz

	while read LINE; do
		REG=$(echo $LINE | awk '{ print $3; }')
		CHUNK_NBR=$(echo $LINE | awk '{ print $1; }')
		OUT=${PREFIX}_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_common.bcf
		LOG=${PREFIX}_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_common.log
		TIM=${PREFIX}_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_common.time
		
		# list of haploid individuals : /mnt/project/data/chrX_males_to_keep_high_het_aneuploid_removed.txt
		dx run app-swiss-army-knife -iimage_file="/${BIN}" --folder="${ODIR}/" -icmd="/usr/bin/time -vo $TIM phase_common_static --input $BCF --map $MAP --output $OUT --thread 72 --log $LOG --filter-maf 0.001 --region $REG --pedigree ${PED} --haploid /mnt/project/data/chrX_males_to_keep_high_het_aneuploid_removed.txt && bcftools index -f $OUT --threads 72" --instance-type mem1_ssd1_v2_x72 --priority normal --name ${PREFIX}_shapeit5_common_chr${CHR}_${CHUNK_NBR} -y --brief --ignore-reuse
		
	done < ${CHUNKS}
done


```
</div>


<br>
#### Ligate chunks
Ligation of the phased common variants is performed using the SHAPEIT5_ligate tool. The program requires an ordered list of the phased common variants. We recommend using appropriate naming of the files in the previous step, so that a command such as ls -1v can directly produce the list of files in the right order.


<div class="code-example" markdown="1">
```bash
#!bin/bash

PREFIX=ukb200k
mkdir -p TMP
BIN=docker/shapeit5.ukb200k.tar.gz
IDIR=Phasing/step1_phase_common/chunks
ODIR=Phasing/step2_ligate
dx mkdir -p ${ODIR}
threads=8

PED=ukb_200k.ped

for CHR in X; do

	SCR=TMP/script_ligate_chr${CHR}.sh
	IN=/mnt/project/${IDIR}/UKB_chr${CHR}.chunk_*.shapeit5_common.bcf
	OUT=${PREFIX}_chr${CHR}.shapeit5_common_ligate.bcf
	printf '#!bin/bash\ndx run app-swiss-army-knife --folder='${ODIR}/' -iimage_file='/${BIN}' ' > ${SCR}
	N=$(($(cat ../chunks/phase_common/chunks_chr${CHR}.txt | wc -l)-1))
	LINE=''

	for SPL in $(seq 0 1 ${N}); do
		LINE="$LINE"-iin="${IDIR}"/"${PREFIX}"_chr"${CHR}".chunk_"${SPL}".shapeit5_common.bcf" "
		LINE="$LINE"-iin="${IDIR}"/"${PREFIX}"_chr"${CHR}".chunk_"${SPL}".shapeit5_common.bcf.csi" "
	done

	LINE2='-iin=data/'${PED}' -icmd="ls -1v '${PREFIX}'_chr'${CHR}'.chunk_*.shapeit5_common.bcf > list_ligate.chr'${CHR}'.txt && ligate_static --pedigree '${PED}' --input list_ligate.chr'${CHR}'.txt --output '${OUT}' --thread '${threads}' --index" --instance-type mem1_ssd1_v2_x8 --priority low --name ligate_chr'${CHR}' -y --brief --ignore-reuse'

	L="$LINE""$LINE2"
	echo ${L} >> ${SCR}
	bash ${SCR}
done
```
</div>

At the end of this step, we have a region-wide haplotype scaffold.

<br>
#### Phasing rare variants in chunks
Having obtained the haplotype scaffold at common variants in the previous step, we can now phase rare variants using the **SHAPEIT5_phase_rare** tool. We always do this in chunks, in order to maximise the parallelization of our jobs (e.g. in 5Mb chunks). The program requires two input files, one is the haplotype scaffold generated in the previous steps, and the other input required is a file containing the whole unphased region, the same input file provided for **SHAPEIT5_phase_common**: SHAPEIT5 automatically discards common variants for this file, duplicated in the phased scaffold.

		
<div class="code-example" markdown="1">
```bash
#!bin/bash


PREFIX=ukb200k
BIN=docker/shapeit5.ukb200k.tar.gz
threads=36
IDIR=Phasing/step2_ligate

ODIR=Phasing/step3_phase_rare/chunks
dx mkdir -p ${ODIR}

PED=/mnt/project/data/ukb_200k.ped


for CHR in X; do
	CHUNKS=../autosomes/chunks/phase_rare/chunks_chr${CHR}.txt
	BCF=/mnt/project/QC/chr${CHR}/ukb24304_chrX_without_PAR.Fmissing.bcf # QCeed data, unphased, full chromosome, full sample set.
	MAP=/mnt/project/data/shapeit_maps/chr${CHR}.b38.gmap.gz
	SCAF=/mnt/project/Phasing/step2_ligate/${PREFIX}_chr${CHR}.shapeit5_common_ligate.bcf

	while read LINE; do
		CHUNK_NBR=$(echo $LINE | awk '{ print $1; }')		
		SCAFFOLD_REG=$(echo $LINE | awk '{ print $3; }')
		SCAFFOLD_REG_START=$(echo ${SCAFFOLD_REG} | cut -d":" -f 2 | cut -d"-" -f1)
		SCAFFOLD_REG_END=$(echo ${SCAFFOLD_REG} | cut -d":" -f 2 | cut -d"-" -f2)
		SCAFFOLD_REG_NAME=${CHR}_${SCAFFOLD_REG_START}_${SCAFFOLD_REG_END}
			
		INPUT_REG=$(echo $LINE | awk '{ print $4; }')
		INPUT_REG_START=$(echo ${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f1)
		INPUT_REG_END=$(echo ${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f2)
		INPUT_REG_NAME=${CHR}_${INPUT_REG_START}_${INPUT_REG_END}

		OUT=${PREFIX}_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_rare.bcf
		LOG=${PREFIX}_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_rare.log
		TIM=${PREFIX}_chr${CHR}.chunk_${CHUNK_NBR}.shapeit5_rare.time

	dx run app-swiss-army-knife -iimage_file="/${BIN}" --folder="${ODIR}/" -icmd="/usr/bin/time -vo $TIM phase_rare_static --input $BCF --map $MAP --output $OUT --thread ${threads} --log $LOG --scaffold $SCAF --scaffold-region $SCAFFOLD_REG --input-region $INPUT_REG --pedigree ${PED} --haploids /mnt/project/data/chrX_males_to_keep_high_het_aneuploid_removed.txt && bcftools index -f $OUT --threads ${threads}" --instance-type mem1_ssd1_v2_x36 --priority normal --name WGS_shapeit5_rare_chr${CHR}_${CHUNK_NBR} -y
		
	done < ${CHUNKS}
done
```
</div>	

The full list of options for the **SHAPEIT5_phase_rare** command can be found [here](https://odelaneau.github.io/shapeit5/docs/documentation/phase_rare/).

<br>

#### Concatenate chunks
Since the phasing of rare variants has been performed using a ligated haplotype scaffold, the chunks can now be concatenated together to resolve the phasing for the entire chromosome. This is done using the **bcftools concat** command.


<div class="code-example" markdown="1">
```bash
#!/bin/bash

PREFIX=ukb200k
threads=16
IDIR=Phasing/step3_phase_rare/chunks
ODIR=Phasing/step3_phase_rare

for CHR in X; do
	chunks=/mnt/project/${IDIR}/${PREFIX}_chr${CHR}.chunk_*.shapeit5_rare.bcf
	OUT=${PREFIX}_chr${CHR}.full.shapeit5_rare.bcf
	dx run app-swiss-army-knife --folder="${ODIR}/" -icmd="ls -1v ${chunks} > concat_list_chr${CHR}.txt && bcftools concat -n -f concat_list_chr${CHR}.txt -o ${OUT} && bcftools index ${OUT} --threads ${threads} && rm concat_list_chr${CHR}.txt" --instance-type mem1_ssd1_v2_x16 --priority normal --name WGS_concat_chr${CHR} -y	
done
```
</div>	



<br>








