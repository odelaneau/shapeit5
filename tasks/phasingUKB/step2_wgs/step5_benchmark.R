
REG=read.table("/home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt", head=FALSE)
TYP=c("fqa","fqc", "fqr")
nTYP=length(TYP)
PAR=c("default", "scaffold", "default.depth8", "scaffold.depth8", "default.modulo5", "scaffold.modulo5")
nPAR=length(PAR)



#BEAGLE
BGLerror = rep(0, nTYP)
BGLtotal = rep(0, nTYP)
for (r in 1:nrow(REG)) {
	for (t in 1:nTYP) {
		tmp=paste("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step3_runbeagle/benchmark_ukb23352_c20_qc_v1.", REG$V3[1], ".beagle5.3.", TYP[t], ".sample.switch.txt.gz", sep="");
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
			tmp=paste("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.", REG$V3[1], ".shapeit5." , PAR[p], ".", TYP[t], ".sample.switch.txt.gz", sep="");
			cat (tmp, "\n")
			SHP = read.table(tmp, head=FALSE)
			SHP = SHP[!is.nan(SHP$V4), ]
			SHPerror[p, t] = SHPerror[p, t] + sum(SHP$V2)
			SHPtotal[p, t] = SHPtotal[p, t] + sum(SHP$V3)
		}
	}
}
	
	

	
#BEAGLE
vBIN=c(0,1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000)
nBIN=length(vBIN)

t=3
A0 = c()
B0 = c()
C0 = c()
for (r in 1:nrow(REG)) {
	tmp=paste("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step3_runbeagle/benchmark_ukb23352_c20_qc_v1.", REG$V3[1], ".beagle5.3.", TYP[t], ".frequency.switch.txt.gz", sep="");
	cat (tmp, "\n")
	BGL = read.table(tmp, head=FALSE)
	if (length(A0) ==0) {
		A0 = BGL$V2
		B0 = BGL$V3
		C0 = BGL$V1
	} else {
		A0 = A0 + BGL$V2
		B0 = B0 + BGL$V3
	}		
}
BIN=cut(C0, breaks=vBIN, labels=1:(nBIN-1))
DF = cbind(A0, B0, BIN)
vE = as.vector(by (DF[, 1], DF[, 3], sum))
vT = as.vector(by (DF[, 2], DF[, 3], sum))
x = vE*100/vT
	
#SHAPEIT
t=3
A1 = c()
B1 = c()
C1 = c()
for (r in 1:nrow(REG)) {
	tmp=paste("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step6_runshapeit/benchmark_ukb23352_c20_qc_v1.", REG$V3[1], ".shapeit5.phase2.bcf.", TYP[t], ".frequency.switch.txt.gz", sep="");
	cat (tmp, "\n")
	BGL = read.table(tmp, head=FALSE)
	if (length(A1) ==0) {
		A1 = BGL$V2
		B1 = BGL$V3
		C1 = BGL$V1
	} else {
		A1 = A1 + BGL$V2
		B1 = B1 + BGL$V3
	}		
}
BIN=cut(C1, breaks=vBIN, labels=1:(nBIN-1))
DF = cbind(A1, B1, BIN)
vE = as.vector(by (DF[, 1], DF[, 3], sum))
vT = as.vector(by (DF[, 2], DF[, 3], sum))
y = vE*100/vT


plot(log10(vBIN[2:nBIN]), x, type="b", col="red", xlab="log10(MAC)", ylab="Switch error (%)")
points(log10(vBIN[2:nBIN]), y, type="b", col="blue")

plot(log10(vBIN[2:nBIN]), (y-x) * 100 / x, type="b", xlab="log10(MAC)", ylab="Switch error reduction (%)", ylim=c(-50, +10))
abline(h=0, col="red")




plot(log10(C0), (A0/B0-A1/B1) / (A0/B0), type="b", xlab="log10(MAC)", ylab="Switch error reduction (%)", ylim=c(-1, +1))
abline(h=0, col="red")






library(RColorBrewer)
COL = brewer.pal(12,"Paired")


pdf ("plot1.pdf", 40,15)

par(mfrow=c(1,2))

D = read.table("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.chr20:7702567-12266861.shapeit5.default.depth8.fqc.block.switch.txt.gz", head=FALSE)
S = read.table("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.chr20:7702567-12266861.shapeit5.default.depth8.fqc.sample.switch.txt.gz", head=FALSE)
S = S[!is.nan(S$V4), ]

SMPL = unique(S$V1)
nSMPL = length(SMPL)
plot(0,0,type="n", xlim=c(min(D$V2), max(D$V2)), ylim=c(0, nSMPL), main="SHAPEIT")

for (i in 1:nSMPL) {
	Ds = D[which(D$V1 == SMPL[i]), ]
	for (l in 2:nrow(Ds)) {
		rect(Ds$V2[l-1], i, Ds$V2[l], i+1, col = ifelse (l %% 2 == 0, COL[1], COL[2]), border = NULL, lwd=0.5)
	}
}


	

D = read.table("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step3_runbeagle/benchmark_ukb23352_c20_qc_v1.chr20:7702567-12266861.beagle5.3.fqc.block.switch.txt.gz", head=FALSE)
S = read.table("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step3_runbeagle/benchmark_ukb23352_c20_qc_v1.chr20:7702567-12266861.beagle5.3.fqc.sample.switch.txt.gz", head=FALSE)
S = S[!is.nan(S$V4), ]

SMPL = unique(S$V1)
nSMPL = length(SMPL)
plot(0,0,type="n", xlim=c(min(D$V2), max(D$V2)), ylim=c(0, nSMPL), main="BEAGLE")

for (i in 1:nSMPL) {
	Ds = D[which(D$V1 == SMPL[i]), ]
	for (l in 2:nrow(Ds)) {
		rect(Ds$V2[l-1], i, Ds$V2[l], i+1, col = ifelse (l %% 2 == 0, COL[3], COL[4]), border = NULL, lwd=0.5)
	}
}


	
dev.off()

