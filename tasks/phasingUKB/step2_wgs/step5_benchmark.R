
REG=read.table("/home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt", head=FALSE)
TYP=c("fqa","fqc")
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
t=2
A = c()
B = c()
C = c()
for (r in 1:nrow(REG)) {
	tmp=paste("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step3_runbeagle/benchmark_ukb23352_c20_qc_v1.", REG$V3[1], ".beagle5.3.", TYP[t], ".frequency.switch.txt.gz", sep="");
	cat (tmp, "\n")
	BGL = read.table(tmp, head=FALSE)
	if (length(A) ==0) {
		A = BGL$V2
		B = BGL$V3
		C = BGL$V1
	} else {
		A = A + BGL$V2
		B = B + BGL$V3
	}		
}
BIN=cut(C, breaks=c(0,1000,2000,5000,10000,20000,50000,100000), labels=1:7)
DF = cbind(A, B, BIN)
vE = as.vector(by (DF[, 1], DF[, 3], sum))
vT = as.vector(by (DF[, 2], DF[, 3], sum))
x = vE*100/vT

	
#SHAPEIT
t=2
A = c()
B = c()
C = c()
for (r in 1:nrow(REG)) {
	tmp=paste("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step4_runshapeit/benchmark_ukb23352_c20_qc_v1.", REG$V3[1], ".shapeit5.", PAR[3], ".", TYP[t], ".frequency.switch.txt.gz", sep="");
	cat (tmp, "\n")
	BGL = read.table(tmp, head=FALSE)
	if (length(A) ==0) {
		A = BGL$V2
		B = BGL$V3
		C = BGL$V1
	} else {
		A = A + BGL$V2
		B = B + BGL$V3
	}		
}
BIN=cut(C, breaks=c(0,1000,2000,5000,10000,20000,50000,100000), labels=1:7)
DF = cbind(A, B, BIN)
vE = as.vector(by (DF[, 1], DF[, 3], sum))
vT = as.vector(by (DF[, 2], DF[, 3], sum))
y = vE*100/vT


plot(log10(c(1000,2000,5000,10000,20000,50000,100000)), x, type="b", col="red")
points(log10(c(1000,2000,5000,10000,20000,50000,100000)), y, type="b", col="blue")


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

