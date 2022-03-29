#!/bin/bash
#REF: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/
for CHR in 20; do
	CHAIN=/mnt/project/Phasing/PhasingSNParray/step4_liftover/hg19ToHg38.over.chain.gz
	VCFF=/mnt/project/Phasing/PhasingSNParray/step3_swapalleles/full_c$CHR\_b0_v2.b37.swapped.vcf.gz
	REF=/mnt/project/Phasing/PhasingSNParray/step4_liftover/GRCh38_full_analysis_set_plus_decoy_hla.fa
	LIFF=full_c$CHR\_b0_v2.b38.vcf.gz
	SORF=full_c$CHR\_b0_v2.b38.sorted.vcf.gz
	dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step4_liftover/" -iimage_file="/docker/liftovervcf_0.0.1.tar.gz" -icmd="liftoverVCF_static --input $VCFF --output $LIFF --chain $CHAIN --fasta $REF --chr chr$CHR && bcftools sort -Oz -m 6G -o $SORF $LIFF && rm $LIFF && bcftools index $SORF" --tag liftover --tag chr$CHR --instance-type mem2_ssd1_v2_x2 --priority normal --name liftover3_chr$CHR -y
done

