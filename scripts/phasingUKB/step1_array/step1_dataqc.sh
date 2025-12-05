#!/bin/bash

#Get SNP list
dx run app-swiss-army-knife -iin="/Bulk/Genotype Results/Genotype calls/ukb_snp_qc.txt" --folder="/Phasing/PhasingSNParray/step1_dataqc/" -icmd="cat ukb_snp_qc.txt | awk '{ print \$1, \$159; }' > SNPlist.unfiltered.txt && cat SNPlist.unfiltered.txt | sed '1d' | awk '{ if (\$2 == 1) print \$1; }' > SNPlist.filtered.QC.txt" --tag qc_snp --instance-type mem1_ssd1_v2_x2 --name qc_snp --priority normal -y

#Get sample lists
dx run app-swiss-army-knife -iin="/Phasing/PhasingSNParray/step1_dataqc/INDlist.unfiltered.txt" -iin="/Phasing/PhasingSNParray/step1_dataqc/samples.parents.txt" -iin="/Phasing/PhasingSNParray/step1_dataqc/samples.related.txt" --folder="/Phasing/PhasingSNParray/step1_dataqc/" -icmd="cat INDlist.unfiltered.txt | sed '1d' | awk '{ if (\$2 == 1) print \$1, \$1; }' > INDlist.filtered.QC.txt && cat INDlist.filtered.QC.txt | grep -vf samples.parents.txt > INDlist.filtered.noparents.txt && cat INDlist.filtered.QC.txt | grep -f samples.related.txt > INDlist.filtered.related.txt" --tag qc_sample --instance-type mem1_ssd1_v2_x2 --name qc_sample --priority normal -y

#Filter each chromosome
for CHR in $(seq 20 20); do
	#FULL DATA
	dx run app-swiss-army-knife -iin="/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c${CHR}_b0_v2.bed" -iin="/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c${CHR}_b0_v2.bim" -iin="/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c${CHR}_b0_v2.fam" -iin="/Phasing/PhasingSNParray/step1_dataqc/INDlist.filtered.QC.txt" -iin="/Phasing/PhasingSNParray/step1_dataqc/SNPlist.filtered.QC.txt" --folder="/Phasing/PhasingSNParray/step1_dataqc/" -icmd="plink2 --bfile ukb22418_c${CHR}_b0_v2 --keep INDlist.filtered.QC.txt --extract SNPlist.filtered.QC.txt --export vcf bgz --out full_c${CHR}_b0_v2 && bcftools index full_c${CHR}_b0_v2.vcf.gz" --tag plink2bcf1 --instance-type mem1_ssd1_v2_x2 --name plink2bcf1 --priority normal -y
done





