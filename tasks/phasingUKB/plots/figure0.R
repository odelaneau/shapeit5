###### COLORS ######

library(RColorBrewer)
COLpair = brewer.pal(12,"Paired")
COLdiff = brewer.pal(8,"Set1")

###### VECTORS ######
SIZE=c(5000,10000,20000,50000,100000,200000,300000,400000,480853)
lSIZE=c("5k","10k","20k","50k","100k","200k","300k","400k","480k")
nSIZE=length(SIZE)
options(scipen=999)

###### DATA ######
D=read.table("../step1_array/data/data.txt", head=TRUE, stringsAsFactors=FALSE)


#########################################################################################
#					FIGURE0						#
#			Plot switch error rates & sepeed vs sample size 		#						
#########################################################################################

pdf("PDFs/figure0.pdf",6,10)
par(mfrow=c(2,1))

plot(0,0, type='n', xlim=c(1, nSIZE), ylim=c(0.0, 1.5), main="A. SNP array - accuracy", xlab="Sample size", ylab="Switch Error rate (%)", xaxt="n")
axis(1, at=1:nSIZE, label=lSIZE, cex.axis=0.8)
abline(h=seq(0.0,1.5,0.1), col="lightgrey", lty=2)
abline(v=1:nSIZE, col="lightgrey", lty=2)
points(1:nSIZE, D$BeagleTrio, type="b", col="black", lwd=2, pch=20)
points(1:nSIZE, D$ShapeitTrio, type="b", col=COLdiff[2], lwd=2, pch=20)
legend("topright", legend=c("Beagle5.4","SHAPEIT5.0"), fill=c("black", COLdiff[2]), bg="white")

plot(0,0, type='n', xlim=c(1, nSIZE), ylim=c(0.0,12), main="B. SNP array - speed", xlab="Sample size", ylab="Wall clock time (h)", xaxt="n")
axis(1, at=1:nSIZE, label=lSIZE, cex.axis=0.8)
abline(h=seq(0,12,1), col="lightgrey", lty=2)
abline(v=1:nSIZE, col="lightgrey", lty=2)
points(1:nSIZE, D$BeagleT/3600, type="b", col="black", lwd=2, pch=20)
points(1:nSIZE, D$ShapeitT/3600, type="b", col=COLdiff[2], lwd=2, pch=20)

dev.off()

