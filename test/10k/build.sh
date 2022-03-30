#!/bin/bash

bcftools view -G -Ob -o msprime.sites.vcf.gz msprime.all.bcf

bcftools view -c 200:minor msprime.sites.vcf.gz | grep -v "#" | awk '{ print $1, $2-1, $2; }' | shuf | head -n 2500 | sort -k2,2n | tr " " "\t" > msprime.array.bed

bcftools view -Ob -o msprime.array.bcf -T msprime.array.bed msprime.all.bcf

bcftools index msprime.array.bcf

rm msprime.sites.vcf.gz msprime.array.bed

