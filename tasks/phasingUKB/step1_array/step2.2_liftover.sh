#!/bin/bash
#REF: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/
for CHR in 20; do
	CHAIN=/mnt/project/Phasing/PhasingSNParray/step2_liftover/hg19ToHg38.over.chain.gz

	VCFB=/mnt/project/Phasing/PhasingSNParray/step2_liftover/benchmark_c$CHR\_b0_v2.b37b.vcf.gz
	VCFF=/mnt/project/Phasing/PhasingSNParray/step2_liftover/full_c$CHR\_b0_v2.b37.vcf.gz
	VCFV=/mnt/project/Phasing/PhasingSNParray/step2_liftover/validation_c$CHR\_b0_v2.b37.vcf.gz
	
	REF=/mnt/project/Phasing/PhasingSNParray/step2_liftover/GRCh38_full_analysis_set_plus_decoy_hla.fa

	OUTB=benchmark_c$CHR\_b0_v2.b38.vcf.gz
	OUTF=full_c$CHR\_b0_v2.b38.vcf.gz
	OUTV=validation_c$CHR\_b0_v2.b38.vcf.gz

	REJB=benchmark_c$CHR\_b0_v2.b38.rejected.vcf.gz
	REJF=full_c$CHR\_b0_v2.b38.rejected.vcf.gz
	REJV=validation_c$CHR\_b0_v2.b38.rejected.vcf.gz
	
	dx run app-swiss-army-knife --folder "/Phasing/PhasingSNParray/step2_liftover/" -icmd="java -jar \$HOME/picard.jar LiftoverVcf -I $VCFV -O $OUTV -C $CHAIN --REJECT $REJV -R $REF --MAX_RECORDS_IN_RAM 1000 --WRITE_ORIGINAL_POSITION --USE_JDK_DEFLATER true --USE_JDK_INFLATER true" --tag liftover --tag chr$CHR --instance-type mem2_ssd1_v2_x32 --priority normal --name liftover_chr$CHR -y
done



