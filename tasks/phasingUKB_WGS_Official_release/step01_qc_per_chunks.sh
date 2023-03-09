#!bin/bash

for CHR in 20; do

	ODIR0=PhasingWGS_Official_release/step0_qc/chunks/support
	dx mkdir -p ${ODIR0}

	JOBID0=$(dx run app-swiss-army-knife -icmd="tabix /mnt/project/Bulk/Whole\ genome\ sequences/Whole\ genome\ GraphTyper\ joint\ call\ pVCF/QC/qc_metrics_graphtyper_v2.7.1_qc.tab.gz chr${CHR} | awk 'BEGIN{print \"##fileformat=VCFv4.2\n##contig=<ID=chr${CHR}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\"}(10^$14<1 && 10^$14>0){AF=10^\$14;MAF=(AF>0.5?1-AF:AF);ExcHet=\$17/(2*(\$16+\$17+\$18)*(AF)*(1-AF));}{if (\$7>0.5 && \$15 <15012 && ExcHet >=0.5 && ExcHet <= 1.5){ split(\$3, array, \":\"); print \$1\"\t\"\$2\"\t\"\$4\"\t\"array[3]\"\t\"array[4]\"\t.\t.\t.\t.\"}}' | bcftools view -Ob -o pass_qc_chr${CHR}.bcf && bcftools index --threads 2 -f pass_qc_chr${CHR}.bcf" --tag filter_tab --tag tabix --tag swiss_army_knife --instance-type mem1_ssd1_v2_x2 --folder="${ODIR0}" --name qc_step03a --priority normal -y | tail -n1 | cut -d" " -f3)

	#the number of chunks is manual but can easily be stored in an array using dx find data for the other chromosomes.
	ODIR1=PhasingWGS_Official_release/step0_qc/chunks
	dx mkdir -p ${ODIR1}
	for CHUNK in {0..1288}; do
		START=$(((CHUNK*50000)+1))
		END=$(((CHUNK+1)*50000))
		JOBID1=$(dx run app-swiss-army-knife -iin="/Bulk/Whole\ genome\ sequences/Whole\ genome\ GraphTyper\ joint\ call\ pVCF/ukb23352_c${CHR}_b${CHUNK}_v1.vcf.gz" -iin="/Bulk/Whole\ genome\ sequences/Whole\ genome\ GraphTyper\ joint\ call\ pVCF/ukb23352_c${CHR}_b${CHUNK}_v1.vcf.gz.tbi" -iin="${ODIR0}/pass_qc_chr${CHR}.bcf" -iin="${ODIR0}/pass_qc_chr${CHR}.bcf.csi" -icmd="bcftools annotate -x ^INFO/AC,^INFO/AN,^FORMAT/GT -Ou ukb23352_c${CHR}_b${CHUNK}_v1.vcf.gz | bcftools view -f PASS -Ou | bcftools norm -m -any -Ou | bcftools view -i 'ALT!=\"*\"' -Ob -o split_ukb23352_c${CHR}_v1.bcf && bcftools index -f split_ukb23352_c${CHR}_v1.bcf && bcftools isec -c none -n=2 -w1 -r chr${CHR}:${START}-${END} split_ukb23352_c${CHR}_v1.bcf pass_qc_chr${CHR}.bcf -Ob -o pass_qc_ukb23352_c${CHR}_${START}_${END}_v1.bcf && bcftools index -f pass_qc_ukb23352_c${CHR}_${START}_${END}_v1.bcf && bcftools +fill-tags -r chr${CHR}:${START}-${END} -Ou pass_qc_ukb23352_c${CHR}_${START}_${END}_v1.bcf -- -t HWE -S /mnt/project/data/ukb_wgs/unphased/qc/support/UKB_samples_with_WGS_british_single_gbr.txt | bcftools view -G -e \"INFO/HWE_GBR < 1e-30\" -Ob -o pass_hwe_chr${CHR}_${START}_${END}.bcf && bcftools index -f pass_hwe_chr${CHR}_${START}_${END}.bcf && bcftools isec -c none -n=2 -w1 -r chr${CHR}:${START}-${END} pass_qc_ukb23352_c${CHR}_${START}_${END}_v1.bcf pass_hwe_chr${CHR}_${START}_${END}.bcf -Ob -o ukb23352_c${CHR}_${START}_${END}_v1.bcf && bcftools index -f ukb23352_c${CHR}_${START}_${END}_v1.bcf && rm -f pass* split_*" --tag qc --tag filter_tab --tag filter_hwe --instance-type mem1_ssd1_v2_x2 --folder="${ODIR1}" --depends-on ${JOBID0} --name qc_chr${CHR}_chunk${CHUNK} --priority normal -y | tail -n1 | cut -d" " -f3)
		
		
	done
done

