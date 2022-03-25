#!/bin/bash
#REF: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/


for CHR in 20; do
	ANN=/mnt/project/Phasing/PhasingSNParray/step2_liftover/chr_rename.txt

	VCFB=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/benchmark_c$CHR\_b0_v2.vcf.gz
	VCFF=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/full_c$CHR\_b0_v2.vcf.gz
	VCFV=/mnt/project/Phasing/PhasingSNParray/step1_dataqc/validation_c$CHR\_b0_v2.vcf.gz

	OUTB=benchmark_c$CHR\_b0_v2.b37.vcf.gz
	OUTF=full_c$CHR\_b0_v2.b37.vcf.gz
	OUTV=validation_c$CHR\_b0_v2.b37.vcf.gz
	
	dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step2_liftover/" -icmd="bcftools annotate -Oz -o $OUTV --rename-chrs $ANN $VCFV && bcftools index $OUTV" --tag updatechr --tag chr$CHR --instance-type mem2_ssd1_v2_x32 --priority normal --name updatechr_chr$CHR -y
done


