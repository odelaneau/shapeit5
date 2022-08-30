###### COLORS ######

library(RColorBrewer)
COLpair = brewer.pal(12,"Paired")
COLdiff = brewer.pal(8,"Set1")

###### VECTORS ######
REG=read.table("/home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt", head=FALSE)
SIZE=c(2000,5000,10000,20000,50000,100000,147754)
BIN=c(0,1,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000, 150000)
PAR=c("default", "scaffold")
TYPE=c("fqa", "fqc", "fqr", "fqf")
nSIZE=length(SIZE)
nTYPE=length(TYPE)
nPAR=length(PAR)
nREG=length(REG)
nBIN=length(BIN)
options(scipen=999)

#########################################################################################
#					FIGURE1A						#
#			Plot switch error rates	vs sample size 				#						
#########################################################################################

#LOAD BEAGLE SER BETWEEN ALL VARIANTS
error = rep(0, 20)
freq = rep(0, 20)
total = rep(0, 20)
for (s in nSIZE:nSIZE) {
	for (r in 1:nrow(REG)) {
		str1="~/Fuse/PhasingUKB/Phasing/PhasingWGS/step5_runshapeit_phase2/N";
		str2="/benchmark_ukb23352_c20_qc_v1.subset.N";
		tmp=paste(str1, SIZE[s], str2, SIZE[s], ".", REG$V3[r], ".shapeit5.default.fqf.calibration.switch.txt.gz", sep="");
		cat (tmp, "\n")
		DATA = read.table(tmp, head=FALSE)
		error = error + DATA$V5
		total = total + DATA$V6
		freq = freq + DATA$V4
		
	}
}

