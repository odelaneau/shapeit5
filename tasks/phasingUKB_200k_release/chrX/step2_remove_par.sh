#!bin/bash



threads=32


ODIR=QC/chrX/PAR_filter/
dx mkdir -p ${ODIR}




# filter QCeed data.
IN=/mnt/project/QC/chrX/ukb24304_chrX.qceed.bcf

REG0=chrX:0-10000
OUT0=ukb24304_chrX_0_10000.bcf
dx run app-swiss-army-knife -icmd="bcftools view -i 'F_MISSING<0.1' --threads ${threads} -r ${REG0} -Ob -o ${OUT0} ${IN} && bcftools index ${OUT0} --threads ${threads}" --instance-type mem2_ssd1_v2_x32 --folder="${ODIR}" --name PAR0_chr${CHR} --priority normal -y --brief



REG1=chrX:2781480-155701384
OUT1=ukb24304_chrX_2781480_155701384.bcf
dx run app-swiss-army-knife -icmd="bcftools view -i 'F_MISSING<0.1' --threads ${threads} -r ${REG1} -Ob -o ${OUT1} ${IN} && bcftools index ${OUT1} --threads ${threads}" --instance-type mem2_ssd1_v2_x32 --folder="${ODIR}" --name PAR1_chr${CHR} --priority normal -y --brief


REG2=chrX:156030896-1000000000
OUT2=ukb24304_chrX_156030896_1000000000.bcf
dx run app-swiss-army-knife -icmd="bcftools view -i 'F_MISSING<0.1' --threads ${threads} -r ${REG2} -Ob -o ${OUT2} ${IN} && bcftools index ${OUT2} --threads ${threads}" --instance-type mem2_ssd1_v2_x32 --folder="${ODIR}" --name PAR2_chr${CHR} --priority normal -y --brief





# merge:
OUT3=ukb24304_chrX_without_PAR.bcf
ODIR3=QC/chrX/
dx run app-swiss-army-knife -icmd="bcftools concat -n -Ob -o ${OUT3} /mnt/project/${ODIR}/${OUT0} /mnt/project/${ODIR}/${OUT1} /mnt/project/${ODIR}/${OUT2}&& bcftools index ${OUT3} --threads ${threads}" --instance-type mem2_ssd1_v2_x32 --folder="${ODIR3}" --name PAR2_chr${CHR} --priority normal -y --brief


OUT4=ukb24304_chrX_without_PAR.sites.bcf
dx run app-swiss-army-knife -icmd="bcftools view -G --threads ${threads} -Ob -o ${OUT4} /mnt/project/${ODIR3}/${OUT3} && bcftools index ${OUT4} --threads ${threads}" --instance-type mem2_ssd1_v2_x32 --folder="${ODIR3}" --name PAR2_chr${CHR} --priority normal -y --brief



