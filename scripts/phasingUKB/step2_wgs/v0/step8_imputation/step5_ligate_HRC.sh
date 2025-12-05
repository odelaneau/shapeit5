#!/bin/bash

#!/bin/bash

METB="hrc"

for CHR in 20; do
	OFOLDER=/Phasing/PhasingWGS/step8_imputation/ligated
	dx run app-swiss-army-knife -iimage_file="docker/shapeit_ligate_0.0.1.tar.gz" -icmd="ls -1v /mnt/project/Phasing/PhasingWGS/step8_imputation/imputation/${METB}/imputed_1k_rp_${METB}_chr${CHR}_*.vcf.gz > list.txt && shapeit_ligate_v0.0.1 --input list.txt --output imputed_1k_rp_${METB}_chr${CHR}.bcf --index --thread 2 && bcftools index -f imputed_1k_rp_${METB}_chr${CHR}.bcf --threads 2 && rm list.txt" --instance-type mem2_ssd1_v2_x2 --priority low --name ${METB}_lg_${CHR} --folder "${OFOLDER}" -y
done

