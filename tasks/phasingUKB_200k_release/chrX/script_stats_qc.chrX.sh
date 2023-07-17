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
