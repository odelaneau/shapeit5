###### COLORS ######
setwd('~/Dropbox/SHAPEIT5/Final_revision/Figures/Scripts/suppFigure2_globalSER/')

library(RColorBrewer)
COLpair = brewer.pal(12,"Paired")
COLdiff = brewer.pal(8,"Set1")

REG=read.table("../../Source_data/chr20.size4Mb.txt", head=FALSE)

SIZE=c(2000,5000,10000,20000,50000,100000,147754)
nSIZE=length(SIZE)
nREG=length(REG)

options(scipen=999)

#LOAD BEAGLE SER BETWEEN ALL VARIANTS / TRIOS+DUOS
# BGLerror_ALL = rep(0, nSIZE)
# BGLtotal_ALL = rep(0, nSIZE)
# for (s in 1:nSIZE) {
# 	for (r in 1:nrow(REG)) {
# 		str1="../../Source_data/WGS/Beagle5.4/N";
# 		str2="/benchmark_ukb23352_c20_qc_v1.subset.N";
# 		str3=".fullchr.";
# 		str4=".fqf.sample.switch.txt.gz"
# 		tmp=paste(str1, SIZE[s], str2, SIZE[s], str3, REG$V4[r], str4, sep="");
# 		cat (tmp, "\n")
# 		BGL = read.table(tmp, head=FALSE)
# 		BGL = BGL[!is.nan(BGL$V4), ]
# 		BGLerror_ALL[s] = BGLerror_ALL[s] + sum(BGL$V2)
# 		BGLtotal_ALL[s] = BGLtotal_ALL[s] + sum(BGL$V3)
# 	}
# }
# 
# write.table(BGLerror_ALL, '../../Source_data/Supp_fig2/BGL_error_ALL.txt', quote=F, col.names=F, row.names=F, sep='\t')
# write.table(BGLtotal_ALL, '../../Source_data/Supp_fig2/BGLtotal_ALL.txt', quote=F, col.names=F, row.names=F, sep='\t')
# 

BGLerror_ALL<-read.table('../../Source_data/Supp_fig2/BGL_error_ALL.txt', hea=F)$V1
BGLtotal_ALL<-read.table('../../Source_data/Supp_fig2/BGLtotal_ALL.txt', hea=F)$V1


# #LOAD BEAGLE SER BETWEEN ALL VARIANTS / TRIOS
# BGLerror_TRIO = rep(0, nSIZE)
# BGLtotal_TRIO = rep(0, nSIZE)
# for (s in 1:nSIZE) {
# 	for (r in 1:nrow(REG)) {
# 		str1="../../Source_data/WGS/Beagle5.4/N";
# 		str2="/benchmark_ukb23352_c20_qc_v1.subset.N";
# 		str3=".fullchr.";
# 		str4=".fqt.sample.switch.txt.gz"
# 		tmp=paste(str1, SIZE[s], str2, SIZE[s], str3, REG$V4[r], str4, sep="");
# 		cat (tmp, "\n")
# 		BGL = read.table(tmp, head=FALSE)
# 		BGL = BGL[!is.nan(BGL$V4), ]
# 		BGLerror_TRIO[s] = BGLerror_TRIO[s] + sum(BGL$V2)
# 		BGLtotal_TRIO[s] = BGLtotal_TRIO[s] + sum(BGL$V3)
# 	}
# }
# write.table(BGLerror_TRIO, '../../Source_data/Supp_fig2/BGLerror_TRIO.txt', quote=F, col.names=F, row.names=F, sep='\t')
# write.table(BGLtotal_TRIO, '../../Source_data/Supp_fig2/BGLtotal_TRIO.txt', quote=F, col.names=F, row.names=F, sep='\t')

BGLerror_TRIO<-read.table('../../Source_data/Supp_fig2/BGLerror_TRIO.txt', hea=F)$V1
BGLtotal_TRIO<-read.table('../../Source_data/Supp_fig2/BGLtotal_TRIO.txt', hea=F)$V1


#LOAD BEAGLE SER BETWEEN COMMON VARIANTS / TRIOS+DUOS
# BGLerror_COM = rep(0, nSIZE)
# BGLtotal_COM = rep(0, nSIZE)
# for (s in 1:nSIZE) {
# 	for (r in 1:nrow(REG)) {
# 		str1="../../Source_data/WGS/Beagle5.4/N";
# 		str2="/benchmark_ukb23352_c20_qc_v1.subset.N";
# 		str3=".fullchr.";
# 		str4=".fqc.sample.switch.txt.gz"
# 		tmp=paste(str1, SIZE[s], str2, SIZE[s], str3, REG$V4[r], str4, sep="");
# 		cat (tmp, "\n")
# 		BGL = read.table(tmp, head=FALSE)
# 		BGL = BGL[!is.nan(BGL$V4), ]
# 		BGLerror_COM[s] = BGLerror_COM[s] + sum(BGL$V2)
# 		BGLtotal_COM[s] = BGLtotal_COM[s] + sum(BGL$V3)
# 	}
# }
# 
# 
# write.table(BGLerror_COM, '../../Source_data/Supp_fig2/BGLerror_COM.txt', quote=F, col.names=F, row.names=F, sep='\t')
# write.table(BGLtotal_COM, '../../Source_data/Supp_fig2/BGLtotal_COM.txt', quote=F, col.names=F, row.names=F, sep='\t')

BGLerror_COM<-read.table('../../Source_data/Supp_fig2/BGLerror_COM.txt', hea=F)$V1
BGLtotal_COM<-read.table('../../Source_data/Supp_fig2/BGLtotal_COM.txt', hea=F)$V1



#LOAD BEAGLE SER BETWEEN ARRAY VARIANTS / TRIOS+DUOS
# BGLerror_ARR = rep(0, nSIZE)
# BGLtotal_ARR = rep(0, nSIZE)
# for (s in 1:nSIZE) {
# 	for (r in 1:nrow(REG)) {
# 		str1="../../Source_data/WGS/Beagle5.4/N";
# 		str2="/benchmark_ukb23352_c20_qc_v1.subset.N";
# 		str3=".fullchr.";
# 		str4=".fqa.sample.switch.txt.gz"
# 		tmp=paste(str1, SIZE[s], str2, SIZE[s], str3, REG$V4[r], str4, sep="");
# 		cat (tmp, "\n")
# 		BGL = read.table(tmp, head=FALSE)
# 		BGL = BGL[!is.nan(BGL$V4), ]
# 		BGLerror_ARR[s] = BGLerror_ARR[s] + sum(BGL$V2)
# 		BGLtotal_ARR[s] = BGLtotal_ARR[s] + sum(BGL$V3)
# 	}
# }
# 
# write.table(BGLerror_ARR, '../../Source_data/Supp_fig2/BGLerror_ARR.txt', quote=F, col.names=F, row.names=F, sep='\t')
# write.table(BGLtotal_ARR, '../../Source_data/Supp_fig2/BGLtotal_ARR.txt', quote=F, col.names=F, row.names=F, sep='\t')

BGLerror_ARR<-read.table('../../Source_data/Supp_fig2/BGLerror_ARR.txt', hea=F)$V1
BGLtotal_ARR<-read.table('../../Source_data/Supp_fig2/BGLtotal_ARR.txt', hea=F)$V1



# 
# 
# #LOAD SHAPEIT SER BETWEEN ALL VARIANTS / TRIOS+DUOS
# SHPerror_ALL = rep(0, nSIZE)
# SHPtotal_ALL = rep(0, nSIZE)
# for (s in 1:nSIZE) {
# 	for (r in 1:nrow(REG)) {
# 		str1="../../Source_data/WGS/Shapeit5/N";
# 		str2="/benchmark_ukb23352_c20_qc_v1.subset.N";
# 		str3=".fullchr.shapeit5.ligated.";
# 		str4=".fqf.sample.switch.txt.gz"
# 		tmp=paste(str1, SIZE[s], str2, SIZE[s], str3, REG$V4[r], str4, sep="");
# 		cat (tmp, "\n")
# 		SHP = read.table(tmp, head=FALSE)
# 		SHP = SHP[!is.nan(SHP$V4), ]
# 		SHPerror_ALL[s] = SHPerror_ALL[s] + sum(SHP$V2)
# 		SHPtotal_ALL[s] = SHPtotal_ALL[s] + sum(SHP$V3)
# 	}
# }
# 
# write.table(SHPerror_ALL, '../../Source_data/Supp_fig2/SHPerror_ALL.txt', quote=F, col.names=F, row.names=F, sep='\t')
# write.table(SHPtotal_ALL, '../../Source_data/Supp_fig2/SHPtotal_ALL.txt', quote=F, col.names=F, row.names=F, sep='\t')

SHPerror_ALL<-read.table('../../Source_data/Supp_fig2/SHPerror_ALL.txt', hea=F)$V1
SHPtotal_ALL<-read.table('../../Source_data/Supp_fig2/SHPtotal_ALL.txt', hea=F)$V1



# 
# 
# #LOAD SHAPEIT SER BETWEEN ALL VARIANTS / TRIOS
# SHPerror_TRIO = rep(0, nSIZE)
# SHPtotal_TRIO = rep(0, nSIZE)
# for (s in 1:nSIZE) {
# 	for (r in 1:nrow(REG)) {
# 		str1="../../Source_data/WGS/Shapeit5/N";
# 		str2="/benchmark_ukb23352_c20_qc_v1.subset.N";
# 		str3=".fullchr.shapeit5.ligated.";
# 		str4=".fqt.sample.switch.txt.gz"
# 		tmp=paste(str1, SIZE[s], str2, SIZE[s], str3, REG$V4[r], str4, sep="");
# 		cat (tmp, "\n")
# 		SHP = read.table(tmp, head=FALSE)
# 		SHP = SHP[!is.nan(SHP$V4), ]
# 		SHPerror_TRIO[s] = SHPerror_TRIO[s] + sum(SHP$V2)
# 		SHPtotal_TRIO[s] = SHPtotal_TRIO[s] + sum(SHP$V3)
# 	}
# }
# 
# write.table(SHPerror_TRIO, '../../Source_data/Supp_fig2/SHPerror_TRIO.txt', quote=F, col.names=F, row.names=F, sep='\t')
# write.table(SHPtotal_TRIO, '../../Source_data/Supp_fig2/SHPtotal_TRIO.txt', quote=F, col.names=F, row.names=F, sep='\t')

SHPerror_TRIO<-read.table('../../Source_data/Supp_fig2/SHPerror_TRIO.txt', hea=F)$V1
SHPtotal_TRIO<-read.table('../../Source_data/Supp_fig2/SHPtotal_TRIO.txt', hea=F)$V1



# 
# #LOAD SHAPEIT SER BETWEEN COMMON VARIANTS / TRIOS+DUOS
# SHPerror_COM = rep(0, nSIZE)
# SHPtotal_COM = rep(0, nSIZE)
# for (s in 1:nSIZE) {
# 	for (r in 1:nrow(REG)) {
# 		str1="../../Source_data/WGS/Shapeit5/N";
# 		str2="/benchmark_ukb23352_c20_qc_v1.subset.N";
# 		str3=".fullchr.shapeit5.ligated.";
# 		str4=".fqc.sample.switch.txt.gz"
# 		tmp=paste(str1, SIZE[s], str2, SIZE[s], str3, REG$V4[r], str4, sep="");
# 		cat (tmp, "\n")
# 		SHP = read.table(tmp, head=FALSE)
# 		SHP = SHP[!is.nan(SHP$V4), ]
# 		SHPerror_COM[s] = SHPerror_COM[s] + sum(SHP$V2)
# 		SHPtotal_COM[s] = SHPtotal_COM[s] + sum(SHP$V3)
# 	}
# }
# 
# 
# write.table(SHPerror_COM, '../../Source_data/Supp_fig2/SHPerror_COM.txt', quote=F, col.names=F, row.names=F, sep='\t')
# write.table(SHPtotal_COM, '../../Source_data/Supp_fig2/SHPtotal_COM.txt', quote=F, col.names=F, row.names=F, sep='\t')

SHPerror_COM<-read.table('../../Source_data/Supp_fig2/SHPerror_COM.txt', hea=F)$V1
SHPtotal_COM<-read.table('../../Source_data/Supp_fig2/SHPtotal_COM.txt', hea=F)$V1



#LOAD SHAPEIT SER BETWEEN ARRAY VARIANTS / TRIOS+DUOS
# SHPerror_ARR = rep(0, nSIZE)
# SHPtotal_ARR = rep(0, nSIZE)
# for (s in 1:nSIZE) {
# 	for (r in 1:nrow(REG)) {
# 		str1="../../Source_data/WGS/Shapeit5/N";
# 		str2="/benchmark_ukb23352_c20_qc_v1.subset.N";
# 		str3=".fullchr.shapeit5.ligated.";
# 		str4=".fqa.sample.switch.txt.gz"
# 		tmp=paste(str1, SIZE[s], str2, SIZE[s], str3, REG$V4[r], str4, sep="");
# 		cat (tmp, "\n")
# 		SHP = read.table(tmp, head=FALSE)
# 		SHP = SHP[!is.nan(SHP$V4), ]
# 		SHPerror_ARR[s] = SHPerror_ARR[s] + sum(SHP$V2)
# 		SHPtotal_ARR[s] = SHPtotal_ARR[s] + sum(SHP$V3)
# 	}
# }
# write.table(SHPerror_ARR, '../../Source_data/Supp_fig2/SHPerror_ARR.txt', quote=F, col.names=F, row.names=F, sep='\t')
# write.table(SHPtotal_ARR, '../../Source_data/Supp_fig2/SHPtotal_ARR.txt', quote=F, col.names=F, row.names=F, sep='\t')

SHPerror_ARR<-read.table('../../Source_data/Supp_fig2/SHPerror_ARR.txt', hea=F)$V1
SHPtotal_ARR<-read.table('../../Source_data/Supp_fig2/SHPtotal_ARR.txt', hea=F)$V1


jpeg("../../Supplementary_Information/Supplementary_figure2.jpeg", 3000, 3000, quality = 100, res=300)
par(mfrow=c(2, 2))
ymax=1.1
plot(0,0, type='n', xlim=c(1, nSIZE), ylim=c(0.0, ymax), main="a. All sites / duos+trios\nWGS N=147,754", xlab="Sample size", ylab="Switch Error rate (%)", xaxt="n")
axis(1, at=1:nSIZE, label=SIZE)
#abline(h=seq(0.0,ymax,0.1), col="lightgrey", lty=3)
#abline(v=1:nSIZE, col="lightgrey", lty=3)
points(1:nSIZE, BGLerror_ALL * 100 / BGLtotal_ALL, type="b", col="black", lwd=2, pch=20)
points(1:nSIZE, SHPerror_ALL * 100 / SHPtotal_ALL, type="b", col=COLdiff[2], lwd=2, pch=20)
legend("topright", fill=c('black',COLdiff[2]), legend=c('Beagle5.4', 'SHAPEIT5.0'), bg="white")


plot(0,0, type='n', xlim=c(1, nSIZE), ylim=c(0.0, ymax), main="b. All sites / trios only\nWGS N=147,754", xlab="Sample size", ylab="Switch Error rate (%)", xaxt="n")
axis(1, at=1:nSIZE, label=SIZE)
#abline(h=seq(0.0,ymax,0.1), col="lightgrey", lty=3)
#abline(v=1:nSIZE, col="lightgrey", lty=3)
points(1:nSIZE, BGLerror_TRIO * 100 / BGLtotal_TRIO, type="b", col="black", lwd=2, pch=20)
points(1:nSIZE, SHPerror_TRIO * 100 / SHPtotal_TRIO, type="b", col=COLdiff[2], lwd=2, pch=20)
legend("topright", fill=c('black',COLdiff[2]), legend=c('Beagle5.4', 'SHAPEIT5.0'), bg="white")


plot(0,0, type='n', xlim=c(1, nSIZE), ylim=c(0.0, 0.5), main="c. Common sites / duos+trios\nWGS N=147,754", xlab="Sample size", ylab="Switch Error rate (%)", xaxt="n")
axis(1, at=1:nSIZE, label=SIZE)
#abline(h=seq(0.0,0.5,0.05), col="lightgrey", lty=3)
#abline(v=1:nSIZE, col="lightgrey", lty=3)
points(1:nSIZE, BGLerror_COM * 100 / BGLtotal_COM, type="b", col="black", lwd=2, pch=20)
points(1:nSIZE, SHPerror_COM * 100 / SHPtotal_COM, type="b", col=COLdiff[2], lwd=2, pch=20)
legend("topright", fill=c('black',COLdiff[2]), legend=c('Beagle5.4', 'SHAPEIT5.0'), bg="white")


plot(0,0, type='n', xlim=c(1, nSIZE), ylim=c(0.0, 2.5), main="d. Axion array sites / duos+trios\nWGS N=147,754", xlab="Sample size", ylab="Switch Error rate (%)", xaxt="n")
axis(1, at=1:nSIZE, label=SIZE)
#abline(h=seq(0.0,2.5,0.5), col="lightgrey", lty=3)
#abline(v=1:nSIZE, col="lightgrey", lty=3)
points(1:nSIZE, BGLerror_ARR * 100 / BGLtotal_ARR, type="b", col="black", lwd=2, pch=20)
points(1:nSIZE, SHPerror_ARR * 100 / SHPtotal_ARR, type="b", col=COLdiff[2], lwd=2, pch=20)
legend("topright", fill=c('black',COLdiff[2]), legend=c('Beagle5.4', 'SHAPEIT5.0'), bg="white")


dev.off()

