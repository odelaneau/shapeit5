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
BGLerror = rep(0, nSIZE)
BGLtotal = rep(0, nSIZE)
for (s in 1:nSIZE) {
	for (r in 1:nrow(REG)) {
		str1="~/Fuse/PhasingUKB/Phasing/PhasingWGS/step3_runbeagle/N";
		str2="/benchmark_ukb23352_c20_qc_v1.subset.N";
		tmp=paste(str1, SIZE[s], str2, SIZE[s], ".", REG$V3[r], ".", TYPE[4], ".sample.switch.txt.gz", sep="");
		cat (tmp, "\n")
		BGL = read.table(tmp, head=FALSE)
		BGL = BGL[!is.nan(BGL$V4), ]
		BGLerror[s] = BGLerror[s] + sum(BGL$V2)
		BGLtotal[s] = BGLtotal[s] + sum(BGL$V3)
	}
}

#LOAD SHAPEIT-DEFAULT SER BETWEEN ALL VARIANTS
SHPerror0 = rep(0, nSIZE)
SHPtotal0 = rep(0, nSIZE)
for (s in 1:nSIZE) {
	for (r in 1:nrow(REG)) {
		str1="~/Fuse/PhasingUKB/Phasing/PhasingWGS/step5_runshapeit_phase2/N";
		str2="/benchmark_ukb23352_c20_qc_v1.subset.N";
		tmp=paste(str1, SIZE[s], str2, SIZE[s], ".", REG$V3[r], ".shapeit5.", PAR[1], ".", TYPE[4], ".sample.switch.txt.gz", sep="");
		cat (tmp, "\n")
		SHP = read.table(tmp, head=FALSE)
		SHP = SHP[!is.nan(SHP$V4), ]
		SHPerror0[s] = SHPerror0[s] + sum(SHP$V2)
		SHPtotal0[s] = SHPtotal0[s] + sum(SHP$V3)
	}
}

#LOAD SHAPEIT-SCAFFOLD SER BETWEEN ALL VARIANTS
SHPerror1 = rep(0, nSIZE)
SHPtotal1 = rep(0, nSIZE)
for (s in 1:nSIZE) {
	for (r in 1:nrow(REG)) {
		str1="~/Fuse/PhasingUKB/Phasing/PhasingWGS/step5_runshapeit_phase2/N";
		str2="/benchmark_ukb23352_c20_qc_v1.subset.N"; 
		tmp=paste(str1, SIZE[s], str2, SIZE[s], ".", REG$V3[r], ".shapeit5.", PAR[2], ".", TYPE[4], ".sample.switch.txt.gz", sep="");
		cat (tmp, "\n")
		SHP = read.table(tmp, head=FALSE)
		SHP = SHP[!is.nan(SHP$V4), ]
		SHPerror1[s] = SHPerror1[s] + sum(SHP$V2)
		SHPtotal1[s] = SHPtotal1[s] + sum(SHP$V3)
	}
}

pdf("PDFs/figure1A.pdf")
plot(0,0, type='n', xlim=c(1, nSIZE), ylim=c(0.35, 1.1), main="All sites", xlab="Sample size", ylab="Switch Error rate (%)", xaxt="n")
axis(1, at=1:nSIZE, label=SIZE)
abline(h=seq(0.3,1.5,0.05), col="lightgrey", lty=2)
abline(v=1:nSIZE, col="lightgrey", lty=2)
points(1:nSIZE, BGLerror * 100 / BGLtotal, type="b", col="black", lwd=2, pch=20)
points(1:nSIZE, SHPerror0 * 100 / SHPtotal0, type="b", col=COLdiff[2], lwd=2, pch=20)
points(1:nSIZE, SHPerror1 * 100 / SHPtotal1, type="b", col=COLdiff[3], lwd=2, pch=20)
legend("topright", legend=c("Beagle5.4","SHAPEIT5.0","SHAPEIT5.0 + scaffold"), fill=c("black", COLdiff[2:3]))
dev.off()


#########################################################################################
#					FIGURE1B						#
#			Plot switch error rates	vs sample size 				#						
#########################################################################################

#LOAD BEAGLE SER BETWEEN ALL VARIANTS
BGLerror = rep(0, nSIZE)
BGLtotal = rep(0, nSIZE)
for (s in 1:nSIZE) {
	for (r in 1:nrow(REG)) {
		str1="~/Fuse/PhasingUKB/Phasing/PhasingWGS/step3_runbeagle/N";
		str2="/benchmark_ukb23352_c20_qc_v1.subset.N";
		tmp=paste(str1, SIZE[s], str2, SIZE[s], ".", REG$V3[r], ".", TYPE[1], ".sample.switch.txt.gz", sep="");
		cat (tmp, "\n")
		BGL = read.table(tmp, head=FALSE)
		BGL = BGL[!is.nan(BGL$V4), ]
		BGLerror[s] = BGLerror[s] + sum(BGL$V2)
		BGLtotal[s] = BGLtotal[s] + sum(BGL$V3)
	}
}

#LOAD SHAPEIT-DEFAULT SER BETWEEN ALL VARIANTS
SHPerror0 = rep(0, nSIZE)
SHPtotal0 = rep(0, nSIZE)
for (s in 1:nSIZE) {
	for (r in 1:nrow(REG)) {
		str1="~/Fuse/PhasingUKB/Phasing/PhasingWGS/step5_runshapeit_phase2/N";
		str2="/benchmark_ukb23352_c20_qc_v1.subset.N";
		tmp=paste(str1, SIZE[s], str2, SIZE[s], ".", REG$V3[r], ".shapeit5.", PAR[1], ".", TYPE[1], ".sample.switch.txt.gz", sep="");
		cat (tmp, "\n")
		SHP = read.table(tmp, head=FALSE)
		SHP = SHP[!is.nan(SHP$V4), ]
		SHPerror0[s] = SHPerror0[s] + sum(SHP$V2)
		SHPtotal0[s] = SHPtotal0[s] + sum(SHP$V3)
	}
}

#LOAD SHAPEIT-SCAFFOLD SER BETWEEN ALL VARIANTS
SHPerror1 = rep(0, nSIZE)
SHPtotal1 = rep(0, nSIZE)
for (s in 1:nSIZE) {
	for (r in 1:nrow(REG)) {
		str1="~/Fuse/PhasingUKB/Phasing/PhasingWGS/step5_runshapeit_phase2/N";
		str2="/benchmark_ukb23352_c20_qc_v1.subset.N"; 
		tmp=paste(str1, SIZE[s], str2, SIZE[s], ".", REG$V3[r], ".shapeit5.", PAR[2], ".", TYPE[1], ".sample.switch.txt.gz", sep="");
		cat (tmp, "\n")
		SHP = read.table(tmp, head=FALSE)
		SHP = SHP[!is.nan(SHP$V4), ]
		SHPerror1[s] = SHPerror1[s] + sum(SHP$V2)
		SHPtotal1[s] = SHPtotal1[s] + sum(SHP$V3)
	}
}

pdf("PDFs/figure1B.pdf")
plot(0,0, type='n', xlim=c(1, nSIZE), ylim=c(0.0, 2.5), main="SNP array sites", xlab="Sample size", ylab="Switch Error rate (%)", xaxt="n")
abline(h=seq(0.0,2.5,0.1), col="lightgrey", lty=2)
abline(v=1:nSIZE, col="lightgrey", lty=2)
axis(1, at=1:nSIZE, label=SIZE)
points(1:nSIZE, BGLerror * 100 / BGLtotal, type="b", col="black", lwd=2, pch=20)
points(1:nSIZE, SHPerror0 * 100 / SHPtotal0, type="b", col=COLdiff[2], lwd=2, pch=20)
points(1:nSIZE, SHPerror1 * 100 / SHPtotal1, type="b", col=COLdiff[3], lwd=2, pch=20)
legend("topright", legend=c("Beagle5.4","SHAPEIT5.0","SHAPEIT5.0 + scaffold"), fill=c("black", COLdiff[2:3]))
dev.off()





