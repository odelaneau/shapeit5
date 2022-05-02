library(RColorBrewer)
COLpair = brewer.pal(12,"Paired")
COLdiff = brewer.pal(8,"Set1")

N=150000
REG=read.table("/home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt", head=FALSE)

#######################################################################################
#
#				PLOT SER BY FREQUENCY
#
#######################################################################################

vBIN=c(0,1,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000, 150000)
nBIN=length(vBIN)

#Computes SER by frequency bin
freqSER <- function (prefix, suffix) {
	A = rep(0, N)
	B = rep(0, N)
	C = 1:N
	for (r in 1:nrow(REG)) {
		tmp=paste(prefix, REG$V3[r], suffix, sep=".");
		DATA = read.table(tmp, head=FALSE)
		cat (tmp, " ", nrow(DATA), " ", max(DATA$V1), "\n")		
		for (l in 1:nrow(DATA)) {
			A[DATA$V1[l]] = A[DATA$V1[l]] + DATA$V2[l]
			B[DATA$V1[l]] = B[DATA$V1[l]] + DATA$V3[l]
		}
	}
	bins=cut(C, breaks=vBIN, labels=1:(nBIN-1))
	merged = cbind(A, B, bins)
	aggrE = as.vector(by(merged[, 1], merged[, 3], sum))
	aggrT = as.vector(by(merged[, 2], merged[, 3], sum))
	ser = aggrE*100/aggrT
	return (ser);
}

#
bgl_default_freq = freqSER("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step3_runbeagle/benchmark_ukb23352_c20_qc_v1", "beagle5.3.fqf.frequency.switch.txt.gz");
shp_default_freq = freqSER("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step6_runshapeit/benchmark_ukb23352_c20_qc_v1", "shapeit5.default.fqf.frequency.switch.txt.gz");
shp_scaffold_freq = freqSER("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step6_runshapeit/benchmark_ukb23352_c20_qc_v1", "shapeit5.scaffold.fqf.frequency.switch.txt.gz");

pdf("figure1.v3.pdf", 12,6)
par(mfrow=c(1,2))
ser_bin=c(0.1, 0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0)

plot(log10(vBIN[2:nBIN]), log10(bgl_default_freq), type="n", pch=20, xlab="MAC", ylab="Switch Error Rate (%)", col="black", lwd=2, xaxt='n', yaxt='n')
abline(h=log10(ser_bin), col="lightgrey", lty=2)
abline(v=log10(vBIN[2:nBIN]), col="lightgrey", lty=2)
axis(1, at=log10(vBIN[2:nBIN]), label=vBIN[2:nBIN], las=2)
axis(2, at=log10(ser_bin), label=ser_bin)
points(log10(vBIN[2:nBIN]), log10(bgl_default_freq), type="o", pch=20, col="black", lwd=2, xaxt='n', yaxt='n')
points(log10(vBIN[2:nBIN]), log10(shp_default_freq), type="o", pch=20, xlab="MAC", col=COLdiff[1], lwd=2)
points(log10(vBIN[2:nBIN]), log10(shp_scaffold_freq), type="o", pch=20, xlab="MAC", col=COLdiff[2], lwd=2)
legend("topright", fill=c("black", COLdiff[1:2]), legend=c("beagle5.3","shapeit5-default","shapeit5-scaffold"), bg="white")
legend("topright", fill=c("black", COLdiff[2]), legend=c("beagle5.3","shapeit5"), bg="white")

plot(log10(vBIN[2:nBIN]), (shp_scaffold_freq - bgl_default_freq) * 100.0 / bgl_default_freq, type="n", pch=20, xlab="MAC", ylab = "Switch Error Rate Reduction (%)", col="black", lwd=2, xaxt='n', ylim=c(-50, 0))
abline(h=seq(-50,0,5), col="lightgrey", lty=2)
abline(v=log10(vBIN[2:nBIN]), col="lightgrey", lty=2)
abline(v=log10(750), col="blue", lty=2)
axis(1, at=log10(vBIN[2:nBIN]), label=vBIN[2:nBIN], las=2)
points(log10(vBIN[2:nBIN]), (shp_default_freq - bgl_default_freq) * 100.0 / bgl_default_freq, type="o", pch=20, col=COLdiff[1], lwd=2)
points(log10(vBIN[2:nBIN]), (shp_scaffold_freq - bgl_default_freq) * 100.0 / bgl_default_freq, type="o", pch=20, col=COLdiff[2], lwd=2)
abline(h=0, col="red")
dev.off()


#######################################################################################
#
#				PLOT SER SCAFFOLD
#
#######################################################################################


TYP=c("fqa","fqc")
nTYP=length(TYP)
PAR=c("default", "scaffold")
nPAR=length(PAR)

#BEAGLE
BGLerror = rep(0, nTYP)
BGLtotal = rep(0, nTYP)
for (r in 1:nrow(REG)) {
	for (t in 1:nTYP) {
		tmp=paste("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step3_runbeagle/benchmark_ukb23352_c20_qc_v1.", REG$V3[r], ".beagle5.3.", TYP[t], ".sample.switch.txt.gz", sep="");
		cat (tmp, "\n")
		BGL = read.table(tmp, head=FALSE)
		BGL = BGL[!is.nan(BGL$V4), ]
		BGLerror[t] = BGLerror[t] + sum(BGL$V2)
		BGLtotal[t] = BGLtotal[t] + sum(BGL$V3)
	}
}
	
#SHAPEIT
SHPerror = matrix(0, nrow=nPAR, ncol=nTYP)
SHPtotal = matrix(0, nrow=nPAR, ncol=nTYP)
for (r in 1:nrow(REG)) {
	for (t in 1:nTYP) {
		for (p in 1:nPAR) {
			tmp=paste("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step6_runshapeit/benchmark_ukb23352_c20_qc_v1.", REG$V3[r], ".shapeit5." , PAR[p], ".bcf.", TYP[t], ".sample.switch.txt.gz", sep="");
			cat (tmp, "\n")
			SHP = read.table(tmp, head=FALSE)
			SHP = SHP[!is.nan(SHP$V4), ]
			SHPerror[p, t] = SHPerror[p, t] + sum(SHP$V2)
			SHPtotal[p, t] = SHPtotal[p, t] + sum(SHP$V3)
		}
	}
}


pdf("figure2.pdf", 4,6)
barplot(cbind(BGLerror * 100.0 / BGLtotal, t(SHPerror* 100.0 / SHPtotal)), ylab="Switch Error Rate (%)", beside=TRUE, names.arg=c("Beagle\ndefault", "Shapeit\ndefault", "Shapeit\nscaffold"), col=COLdiff[1:2])
legend("topright", legend=c("Assay sites", "All sites"), fill=COLdiff)
dev.off()

#######################################################################################
#
#				PLOT SWITCHES ARRAY
#
#######################################################################################

library(RColorBrewer)
COL = brewer.pal(12,"Paired")


pdf ("figure3.pdf", 30,30)

par(mfrow=c(1,3), mar=c(4.1,0,4.1,0))


D = read.table("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step3_runbeagle/benchmark_ukb23352_c20_qc_v1.chr20:7702567-12266861.beagle5.3.fqa.block.switch.txt.gz", head=FALSE)
S = read.table("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step3_runbeagle/benchmark_ukb23352_c20_qc_v1.chr20:7702567-12266861.beagle5.3.fqa.sample.switch.txt.gz", head=FALSE)
S = S[!is.nan(S$V4), ]

SMPL = unique(S$V1)
nSMPL = length(SMPL)
plot(0,0,type="n", xlim=c(min(D$V2), max(D$V2)), ylim=c(0, nSMPL), main="Beagle\ndefault", yaxt='n')

for (i in 1:nSMPL) {
	Ds = D[which(D$V1 == SMPL[i]), ]
	for (l in 2:nrow(Ds)) {
		rect(Ds$V2[l-1], i, Ds$V2[l], i+1, col = ifelse (l %% 2 == 0, COL[1], COL[2]), border = NULL, lwd=0.5)
	}
}



D = read.table("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step6_runshapeit/benchmark_ukb23352_c20_qc_v1.chr20:7702567-12266861.shapeit5.default.bcf.fqa.block.switch.txt.gz", head=FALSE)
S = read.table("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step6_runshapeit/benchmark_ukb23352_c20_qc_v1.chr20:7702567-12266861.shapeit5.default.bcf.fqa.sample.switch.txt.gz", head=FALSE)
S = S[!is.nan(S$V4), ]

SMPL = unique(S$V1)
nSMPL = length(SMPL)
plot(0,0,type="n", xlim=c(min(D$V2), max(D$V2)), ylim=c(0, nSMPL), main="Shapeit\ndefault", yaxt='n')

for (i in 1:nSMPL) {
	Ds = D[which(D$V1 == SMPL[i]), ]
	for (l in 2:nrow(Ds)) {
		rect(Ds$V2[l-1], i, Ds$V2[l], i+1, col = ifelse (l %% 2 == 0, COL[3], COL[4]), border = NULL, lwd=0.5)
	}
}


D = read.table("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step6_runshapeit/benchmark_ukb23352_c20_qc_v1.chr20:7702567-12266861.shapeit5.scaffold.bcf.fqa.block.switch.txt.gz", head=FALSE)
S = read.table("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step6_runshapeit/benchmark_ukb23352_c20_qc_v1.chr20:7702567-12266861.shapeit5.scaffold.bcf.fqa.sample.switch.txt.gz", head=FALSE)
S = S[!is.nan(S$V4), ]

SMPL = unique(S$V1)
nSMPL = length(SMPL)
plot(0,0,type="n", xlim=c(min(D$V2), max(D$V2)), ylim=c(0, nSMPL), main="Shapeit\nscaffold", yaxt='n')

for (i in 1:nSMPL) {
	Ds = D[which(D$V1 == SMPL[i]), ]
	for (l in 2:nrow(Ds)) {
		rect(Ds$V2[l-1], i, Ds$V2[l], i+1, col = ifelse (l %% 2 == 0, COL[5], COL[6]), border = NULL, lwd=0.5)
	}
}

dev.off()

