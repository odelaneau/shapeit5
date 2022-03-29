#!/bin/bash

BCF=/mnt/project/data/ukb_wgs/unphased/qc/ukb23352_c20_qc_v1.bcf

#GET VALIDATION DATA
SUB=/mnt/project/Phasing/PhasingWGS/step1_preparedata/validation.samples.txt
OUT=validation_ukb23352_c20_qc_v1.bcf
dx run app-swiss-army-knife --folder "/Phasing/PhasingWGS/step1_preparedata/" -icmd="bcftools view -Ob -o $OUT -S $SUB --force-samples $BCF && bcftools index $OUT" --tag subsetV --tag benchWGS --instance-type mem2_ssd1_v2_x2 --name benchWGS_subsetV --priority normal -y

#GET BENCHMARK DATA
SUB=/mnt/project/Phasing/PhasingWGS/step1_preparedata/INDlist.filtered.noparents.txt
OUT=benchmark_ukb23352_c20_qc_v1.bcf
dx run app-swiss-army-knife --folder "/Phasing/PhasingWGS/step1_preparedata/" -icmd="bcftools view -Ob -o $OUT -S $SUB --force-samples $BCF && bcftools index $OUT" --tag subsetB --tag benchWGS --instance-type mem2_ssd1_v2_x2 --name benchWGS_subsetB --priority normal -y

