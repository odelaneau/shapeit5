setwd('~/Dropbox/SHAPEIT5/Final_revision/Figures/Scripts/extended_figure1/')
#pdf("Extended_figure1.pdf", 10, 5)
jpeg("../../Extended_figures/Extended_figure1.jpeg", 3500, 1700, quality = 100, res=300)

par(mfrow=c(1,2))

#########################################################################################
#				FIGURE 2A and 2B					#
#			SER and RT vs sample size SNP ARRAY				#
#########################################################################################

library(RColorBrewer)
COLpair = brewer.pal(12,"Paired")
COLdiff = brewer.pal(8,"Set1")
SIZE=c(5000,10000,20000,50000,100000,200000,300000,400000,480853)
lSIZE=c("5k","10k","20k","50k","100k","200k","300k","400k","480k")
nSIZE=length(SIZE)
options(scipen=999)
D=read.table("../../Source_data/Extended_fig1/", head=TRUE, stringsAsFactors=FALSE)



plot(0,0, type='n', xlim=c(1, nSIZE), ylim=c(0.0, 1.5), main="a. SNP ARRAY - accuracy", xlab="Sample size", ylab="Switch Error rate (%)", xaxt="n")
text(1:nSIZE, par("usr")[3], labels = lSIZE, srt = 45, adj = c(1.1,1.1), xpd = TRUE) 
#abline(h=seq(0.0,1.5,0.1), col="lightgrey", lty=3)
#abline(v=1:nSIZE, col="lightgrey", lty=3)
points(1:nSIZE, D$BeagleTrio, type="b", col="black", lwd=2, pch=20)
points(1:nSIZE, D$ShapeitTrio, type="b", col=COLdiff[2], lwd=2, pch=20)

legend("topright", legend=c("Beagle5.4","SHAPEIT5.0"), fill=c("black", COLdiff[2]), bg="white")

plot(0,0, type='n', xlim=c(1, nSIZE), ylim=c(0.0,12), main="b. SNP ARRAY - timing", xlab="Sample size", ylab="Wall clock time (h)", xaxt="n")
text(1:nSIZE, par("usr")[3], labels = lSIZE, srt = 45, adj = c(1.1,1.1), xpd = TRUE) 
#abline(h=seq(0,12,1), col="lightgrey", lty=3)
#abline(v=1:nSIZE, col="lightgrey", lty=3)
points(1:nSIZE, D$BeagleT/3600, type="b", col="black", lwd=2, pch=20)
points(1:nSIZE, D$ShapeitT/3600, type="b", col=COLdiff[2], lwd=2, pch=20)
legend("topleft", legend=c("Beagle5.4","SHAPEIT5.0"), fill=c("black", COLdiff[2]), bg="white")


dev.off()

