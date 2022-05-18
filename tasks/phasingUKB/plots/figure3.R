library(RColorBrewer)
COL = brewer.pal(12,"Paired")



D0 = read.table("~/Fuse/PhasingUKB/Phasing/PhasingSNParray/step5_benchmark/benchmark_c20_b0_v2.b38.sorted.N480853.phased.fqa.block.switch.txt.gz", head=FALSE)
D1 = read.table("~/Fuse/PhasingUKB/Phasing/PhasingSNParray/step5_benchmark/benchmark_c20_b0_v2.b38.sorted.N480853.phased.wgs.fqa.block.switch.txt.gz", head=FALSE)
S = read.table("~/Fuse/PhasingUKB/Phasing/PhasingSNParray/step5_benchmark/benchmark_c20_b0_v2.b38.sorted.N480853.phased.wgs.fqa.sample.mendel.txt.gz", head=FALSE, stringsAsFactors=FALSE)
S = S[which(S$V2>=0 & S$V3>=0), ]

SMPL = unique(S$V1)
nSMPL = length(SMPL)

D0 = D0[which(D0$V1 %in% SMPL), ]
D1 = D1[which(D1$V1 %in% SMPL), ]

pdf ("PDFs/figure3.pdf", 9, 9)

par(mfrow=c(2,1), mar=c(4.1,1.1,2.1,1.1))
plot(0,0,type="n", xlim=c(min(D0$V2), max(D0$V2)), ylim=c(0, nSMPL+1), main="SHAPEIT5 - Array validation", yaxt='n', xlab="", ylab="")
for (i in 1:nSMPL) {
	Ds = D0[which(D0$V1 == SMPL[i]), ]
	for (l in 2:nrow(Ds)) {
		rect(Ds$V2[l-1], i, Ds$V2[l], i+1, col = ifelse (l %% 2 == 0, COL[1], COL[2]), border = NULL, lwd=0.5)
	}
}

plot(0,0,type="n", xlim=c(min(D1$V2), max(D1$V2)), ylim=c(0, nSMPL+1), main="SHAPEIT5 - WGS validation", yaxt='n', xlab="Genomic coordinates on chr20", ylab="")
for (i in 1:nSMPL) {
	Ds = D1[which(D1$V1 == SMPL[i]), ]
	for (l in 2:nrow(Ds)) {
		rect(Ds$V2[l-1], i, Ds$V2[l], i+1, col = ifelse (l %% 2 == 0, COL[1], COL[2]), border = NULL, lwd=0.5)
	}
}

dev.off()


