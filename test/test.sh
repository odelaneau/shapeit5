#!/bash/bin

#step1: phasing common variants
./phase_common/bin/SHAPEIT5_phase_common_static --input 10k/msprime.nodup.bcf --filter-maf 0.001  --output 10k/msprime.common.phased.bcf --region 1 --thread 8
bcftools index 10k/msprime.common.phased.bcf

#TODO: make the 10 Mb region split into two overlapping regions of 6Mb

#TODO: phase each chunk independently

#TODO: ligate two chunks together

#step2: validation of haplotypes at common variants
../switch/bin/SHAPEIT5_switch_static --validation 10k/msprime.nodup.bcf --estimation 10k/msprime.common.phased.bcf --region 1 --output 10k/msprime.common.phased

#step3: example of how to phase rare variants ub a 1Mb region
./phase_rare/bin/SHAPEIT5_phase_rare_static --input-plain 10k/msprime.nodup.bcf --scaffold 10k/msprime.common.truth.bcf --output 10k/msprime.rare.chunk1.bcf --scaffold-region 1:1000000-3000000 --input-region 1:1500000-2500000 --thread 8
bcftools index 10k/msprime.rare.test1.bcf

#TODO: bcftools concat all chunks together

#step4: validate rare variants
../switch/bin/SHAPEIT5_switch_static --validation 10k/msprime.nodup.bcf --estimation 10k/msprime.rare.chunk1.bcf --region 1 --output 10k/msprime.rare.chunk1












