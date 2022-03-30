#!/bash/bin

#step1: phasing array data
../phase1/bin/SHAPEIT5_phase1 --input 10k/msprime.array.bcf --output 10k/msprime.array.phased.bcf --region 1 --thread 8
bcftools index 10k/msprime.array.phased.bcf

#step2: validation of array haplotypes
../switch/bin/SHAPEIT5_switch --validation 10k/msprime.all.bcf --estimation 10k/msprime.array.phased.bcf --region 1 --output 10k/msprime.array.phase

#step3: phasing common variants onto array based scafold
../phase1/bin/SHAPEIT5_phase1 --input 10k/msprime.all.bcf --scaffold 10k/msprime.array.phased.bcf --filter-maf 0.001  --output 10k/msprime.scaffold.phased.bcf --region 1 --thread 8
bcftools index 10k/msprime.scaffold.phased.bcf

../phase1/bin/SHAPEIT5_phase1 --input 10k/msprime.all.bcf --filter-maf 0.001  --output 10k/msprime.scaffold.phased.bcf --region 1 --thread 8
bcftools index 10k/msprime.scaffold.phased.bcf


#step4: validation of haplotypes at common variants
../switch/bin/SHAPEIT5_switch --validation 10k/msprime.all.bcf --estimation 10k/msprime.scaffold.phased.bcf --region 1 --output 10k/msprime.scaffold.phased









