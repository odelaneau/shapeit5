#!/bin/bash

#Setting 1
PATHB="Phasing/PhasingWGS/step3_runbeagle/N147754"
EXTB="vcf.gz"
METB="beagle5.4"

#Setting 2
PATHS="Phasing/PhasingWGS/step5_runshapeit_phase2/N147754"
EXTS="default.bcf"
METS="shapeit5"

for CHR in 20; do
		#Setting 2PhasingWGS_Official_release/step2_phase_rare/concat_benchmark/
		#JOBID0=$(dx run app-swiss-army-knife -iin="/Phasing/PhasingWGS/step8_imputation/support/1k_samples_wbi.txt" -icmd="bcftools view /mnt/project/PhasingWGS_Official_release/step3_phase_beagle/UKB_chr${CHR}.full.${METB}.vcf.gz -r chr${CHR} -S ^1k_samples_wbi.txt -Ou  --threads 4 | bcftools filter -e 'AC==0 || AC==AN' -Oz -o UKB_chr${CHR}.rp.${METB}.vcf.gz --threads 4 && bcftools index -f UKB_chr${CHR}.rp.${METB}.vcf.gz --threads 4" --folder="/PhasingWGS_Official_release/step3_phase_beagle/reference_panels" --instance-type mem1_ssd1_v2_x4 --priority low --name subset_${METB}_chr${CHR} -y | tail -n1 | cut -d" " -f3) 

		#JOBID1=$(dx run app-swiss-army-knife -iin="/docker/bref3.22Jul22.46e.jar" -icmd="java -Xmx32G -jar bref3.22Jul22.46e.jar /mnt/project/PhasingWGS_Official_release/step3_phase_beagle/reference_panels/UKB_chr${CHR}.rp.${METB}.vcf.gz > UKB_chr${CHR}.rp.${METB}.bref3" --folder="/PhasingWGS_Official_release/step3_phase_beagle/reference_panels" --instance-type mem3_ssd1_v2_x4 --priority high --name bref_${METB}_chr${CHR} -y | tail -n1 | cut -d" " -f3) #--depends-on ${JOBID0}

		#from vcf.gz to bref3	
		#JOBID2=$(dx run app-swiss-army-knife -iin="/docker/beagle.22Jul22.46e.jar" -icmd="java -Xmx32G -jar beagle.22Jul22.46e.jar gt=/mnt/project/Phasing/PhasingWGS/step8_imputation/target_snp_array/axiom_1k_c${CHR}_b0_v2.b38.vcf.gz map=/mnt/project/data/plink_maps/plink.prefix.chr${CHR}.GRCh38.map ref=/mnt/project/PhasingWGS_Official_release/step3_phase_beagle/reference_panels/UKB_chr${CHR}.rp.${METB}.bref3 out=imputed_1k_rp_${METB}_chr${CHR} nthreads=4 chrom=chr${CHR} impute=true gp=true && bcftools index -f imputed_1k_rp_${METB}_chr${CHR}.vcf.gz --threads 4" --folder="/PhasingWGS_Official_release/step3_phase_beagle/imputation" --instance-type mem3_ssd1_v2_x4 --priority low --name imp_rp_${METB}_chr${CHR} -y --depends-on ${JOBID1} | tail -n1 | cut -d" " -f3)

		####################################

		#Setting 2
		JOBID0=$(dx run app-swiss-army-knife -iin="/Phasing/PhasingWGS/step8_imputation/support/1k_samples_wbi.txt" -icmd="bcftools view /mnt/project/PhasingWGS_Official_release/step2_phase_rare/concat_benchmark/chr${CHR}/UKB_chr${CHR}.full.${METS}_rare.bcf -r chr${CHR} -S ^1k_samples_wbi.txt -Ou  --threads 4 | bcftools filter -e 'AC==0 || AC==AN' -Oz -o UKB_chr${CHR}.rp.${METS}.vcf.gz --threads 4 && bcftools index -f UKB_chr${CHR}.rp.${METS}.vcf.gz --threads 4" --folder="/PhasingWGS_Official_release/step2_phase_rare/reference_panels" --instance-type mem1_ssd1_v2_x4 --priority low --name subset_${METS}_chr${CHR} -y | tail -n1 | cut -d" " -f3) 

		JOBID1=$(dx run app-swiss-army-knife -iin="/docker/bref3.22Jul22.46e.jar" -icmd="java -Xmx32G -jar bref3.22Jul22.46e.jar /mnt/project/PhasingWGS_Official_release/step2_phase_rare/reference_panels/UKB_chr${CHR}.rp.${METS}.vcf.gz > UKB_chr${CHR}.rp.${METS}.bref3" --folder="/PhasingWGS_Official_release/step2_phase_rare/reference_panels" --instance-type mem3_ssd1_v2_x4 --depends-on ${JOBID0} --priority high --name bref_${METS}_chr${CHR} -y | tail -n1 | cut -d" " -f3) #--depends-on ${JOBID0}

		#from vcf.gz to bref3	
		JOBID2=$(dx run app-swiss-army-knife -iin="/docker/beagle.22Jul22.46e.jar" -icmd="java -Xmx32G -jar beagle.22Jul22.46e.jar gt=/mnt/project/Phasing/PhasingWGS/step8_imputation/target_snp_array/axiom_1k_c${CHR}_b0_v2.b38.vcf.gz map=/mnt/project/data/plink_maps/plink.prefix.chr${CHR}.GRCh38.map ref=/mnt/project/PhasingWGS_Official_release/step2_phase_rare/reference_panels/UKB_chr${CHR}.rp.${METS}.bref3 out=imputed_1k_rp_${METS}_chr${CHR} nthreads=4 chrom=chr${CHR} impute=true gp=true && bcftools index -f imputed_1k_rp_${METS}_chr${CHR}.vcf.gz --threads 4" --folder="/PhasingWGS_Official_release/step2_phase_rare/imputation" --instance-type mem3_ssd1_v2_x4 --priority low --name imp_rp_${METS}_chr${CHR} -y --depends-on ${JOBID1} | tail -n1 | cut -d" " -f3)

done
