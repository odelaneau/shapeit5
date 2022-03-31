#!/bash/bin

#step1: phasing array data
../phase1/bin/SHAPEIT5_phase1 --input 10k/msprime.array.bcf --output 10k/msprime.array.phased.bcf --region 1 --thread 8
bcftools index 10k/msprime.array.phased.bcf

#step2: validation of array haplotypes
../switch/bin/SHAPEIT5_switch --validation 10k/msprime.nodup.bcf --estimation 10k/msprime.array.phased.bcf --region 1 --output 10k/msprime.array.phase

#step3.0: phasing common variants 
../phase1/bin/SHAPEIT5_phase1 --input 10k/msprime.nodup.bcf --filter-maf 0.001  --output 10k/msprime.nscaffold.phased.bcf --region 1 --thread 8
bcftools index 10k/msprime.nscaffold.phased.bcf

#step3.1: phasing common variants onto array based scaffold
../phase1/bin/SHAPEIT5_phase1 --input 10k/msprime.nodup.bcf --scaffold 10k/msprime.array.phased.bcf --filter-maf 0.001 --output 10k/msprime.yscaffold.phased.bcf --region 1 --thread 8
bcftools index 10k/msprime.yscaffold.phased.bcf

#step4: validation of haplotypes at common variants
../switch/bin/SHAPEIT5_switch --validation 10k/msprime.nodup.bcf --estimation 10k/msprime.nscaffold.phased.bcf --region 1 --output 10k/msprime.nscaffold.phased
../switch/bin/SHAPEIT5_switch --validation 10k/msprime.nodup.bcf --estimation 10k/msprime.yscaffold.phased.bcf --region 1 --output 10k/msprime.yscaffold.phased

#step5: phase rare variants
../phase2/bin/SHAPEIT5_phase2 --input 10k/msprime.nodup.bcf --scaffold 10k/msprime.nscaffold.phased.bcf --output 10k/msprime.nodup.nscaffold.phased.bcf --scaffold-region 1:1000000-3000000 --input-region 1:1500000-2500000 --thread 8
../phase2/bin/SHAPEIT5_phase2 --input 10k/msprime.nodup.bcf --scaffold 10k/msprime.yscaffold.phased.bcf --output 10k/msprime.nodup.yscaffold.phased.bcf --scaffold-region 1:1000000-3000000 --input-region 1:1500000-2500000 --thread 8
bcftools index 10k/msprime.nodup.nscaffold.phased.bcf
bcftools index 10k/msprime.nodup.yscaffold.phased.bcf

#step6: validate rare variants
../switch/bin/SHAPEIT5_switch --validation 10k/msprime.nodup.bcf --estimation 10k/msprime.nodup.nscaffold.phased.bcf --region 1 --output 10k/msprime.nodup.nscaffold.phased
../switch/bin/SHAPEIT5_switch --validation 10k/msprime.nodup.bcf --estimation 10k/msprime.nodup.yscaffold.phased.bcf --region 1 --output 10k/msprime.nodup.yscaffold.phased

#step7: cleaning
rm 10k/*.gz
rm 10k/*phased*












