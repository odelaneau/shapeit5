#!/bin/bash

bcftools view -G msprime.all.bcf | grep -v "#" | cut -f1-5 | sort | uniq -c | awk '{ if ($1 > 1) print $2, $3-1, $3; }' | sort -k2,2n | tr " " "\t" > duplicates.bed

bcftools view -Ob -o msprime.nodup.bcf -T ^duplicates.bed msprime.all.bcf

bcftools index msprime.nodup.bcf

bcftools view -G -Ob -o msprime.sites.vcf.gz msprime.nodup.bcf

#bcftools view -c 20:minor msprime.sites.vcf.gz | grep -v "#" | awk '{ print $1, $2-1, $2; }' | tr " " "\t" > msprime.common.bed
bcftools view -c 200:minor msprime.sites.vcf.gz | grep -v "#" | awk '{ print $1, $2-1, $2; }' | tr " " "\t" > msprime.common.bed

bcftools view -Ob -o msprime.common.truth.bcf -T msprime.common.bed msprime.nodup.bcf

bcftools index msprime.common.truth.bcf

rm msprime.sites.vcf.gz msprime.common.bed duplicates.bed

