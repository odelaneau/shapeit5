#!/bin/bash


dx run app-swiss-army-knife -iin="/Bulk/Genotype Results/Genotype calls/ukb_snp_qc.txt" --folder="/PhasingSNParray/" -icmd="cat ukb_snp_qc.txt | awk '{ print \$1; }' > test2.txt" --tag qc1 --instance-type mem1_ssd1_v2_x2 --name qc_1 --priority normal -y 

