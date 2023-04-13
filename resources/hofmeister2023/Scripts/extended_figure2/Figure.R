setwd('~/Dropbox/SHAPEIT5/Final_revision/Figures/Scripts/extended_figure2/')

library(RColorBrewer)
COL = brewer.pal(12,"Paired")

D0<-read.table('../../Source_data/SER/S5_array.txt', hea=F)
D1<-read.table('../../Source_data/SER/S5_wgs.txt', hea=F)
SMPL = unique(D0$V1)
nSMPL = length(SMPL)


#pdf ("extended_figure2.pdf", 9, 9)
jpeg("../../Extended_figures/Extended_figure2.jpeg", 2800, 3000, quality = 100, res=300)

par(mfrow=c(2,1), mar=c(4.1,1.1,2.1,1.1))
plot(0,0,type="n", xlim=c(min(D0$V2), max(D0$V2)), ylim=c(0, nSMPL+1), main="SHAPEIT5 - Array validation", yaxt='n', xlab="", ylab="", xaxt='n')
axis(1, at = c(0,1E7,2E7,3E7,4E7,5E7,6E7 ), labels = c("0e+00", "1e+07", "2e+07", "3e+07","4e+07","5e+07","6e+07"))
for (i in 1:nSMPL) {
	Ds = D0[which(D0$V1 == SMPL[i]), ]
	for (l in 2:nrow(Ds)) {
		rect(Ds$V2[l-1], i, Ds$V2[l], i+1, col = ifelse (l %% 2 == 0, COL[1], COL[2]), border = NULL, lwd=0.5)
	}
}

plot(0,0,type="n", xlim=c(min(D1$V2), max(D1$V2)), ylim=c(0, nSMPL+1), main="SHAPEIT5 - WGS validation", yaxt='n', xlab="Genomic coordinates on chr20", ylab="", xaxt='n')
axis(1, at = c(0,1E7,2E7,3E7,4E7,5E7,6E7 ), labels = c("0e+00", "1e+07", "2e+07", "3e+07","4e+07","5e+07","6e+07"))
for (i in 1:nSMPL) {
	Ds = D1[which(D1$V1 == SMPL[i]), ]
	for (l in 2:nrow(Ds)) {
		rect(Ds$V2[l-1], i, Ds$V2[l], i+1, col = ifelse (l %% 2 == 0, COL[1], COL[2]), border = NULL, lwd=0.5)
	}
}

dev.off()


