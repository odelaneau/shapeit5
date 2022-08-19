#!/bash/bin

#step1: phasing common variants 
#../phase1/bin/SHAPEIT5_phase1 --input 10k/msprime.nodup.bcf --filter-maf 0.001  --output 10k/msprime.common.phased.bcf --region 1 --thread 8
#bcftools index 10k/msprime.common.phased.bcf

#step2: validation of haplotypes at common variants
#../switch/bin/SHAPEIT5_switch --validation 10k/msprime.nodup.bcf --estimation 10k/msprime.common.phased.bcf --region 1 --output 10k/msprime.common.phased

#step3: phase rare variants
../phase2/bin/SHAPEIT5_phase2 --input 10k/msprime.nodup.bcf --scaffold 10k/msprime.common.truth.bcf --output 10k/msprime.rare.test1.bcf --scaffold-region 1:1000000-3000000 --input-region 1:1500000-2500000 --thread 8
bcftools index 10k/msprime.rare.test1.bcf

#step4: validate rare variants
../switch/bin/SHAPEIT5_switch --validation 10k/msprime.nodup.bcf --estimation 10k/msprime.rare.test1.bcf --region 1 --output 10k/msprime.rare.test1

#step7: cleaning
rm 10k/*.gz
rm 10k/*phased*












