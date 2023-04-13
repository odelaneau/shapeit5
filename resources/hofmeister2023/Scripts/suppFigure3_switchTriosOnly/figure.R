setwd('~/Dropbox/SHAPEIT5/Final_revision/Figures/Scripts/suppFigure3_switchTriosOnly/')

###### COLORS ######
library(RColorBrewer)
COLpair = brewer.pal(12,"Paired")
COLdiff = brewer.pal(8,"Set1")

put.fig.letter <- function(label, location="topleft", x=NULL, y=NULL, 
                           offset=c(0, 0), ...) {
  if(length(label) > 1) {
    warning("length(label) > 1, using label[1]")
  }
  if(is.null(x) | is.null(y)) {
    coords <- switch(location,
                     topleft = c(0.025,0.97),
                     topcenter = c(0.5525,0.98),
                     topright = c(0.985, 0.98),
                     bottomleft = c(0.025, 0.02), 
                     bottomcenter = c(0.5525, 0.02), 
                     bottomright = c(0.985, 0.02),
                     c(0.015, 0.98) )
  } else {
    coords <- c(x,y)
  }
  this.x <- grconvertX(coords[1] + offset[1], from="nfc", to="user")
  this.y <- grconvertY(coords[2] + offset[2], from="nfc", to="user")
  text(labels=label[1], x=this.x, y=this.y, xpd=T, ...)
}



#########################################################################################
#		LOAD DATA FOR FIGURE 2A							#
#########################################################################################

REG=read.table("../../Source_data/chr20.size4Mb.txt", head=FALSE)
BIN=c(0,1,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000, 100000000)
lBIN=c("singleton","2-5","6-10","11-20","21-50","51-100","101-200","201-500","501-1k","1k-2k","2k-5k","5k-10k","10k-20k","20k-50k","50k-100k", "100k+")
nREG=length(REG)
nBIN=length(BIN)
options(scipen=999)

freqSER0 <- function (prefix, suffix, n) {
  A = rep(0, n)
  B = rep(0, n)
  C = 1:n
  for (r in 1:nrow(REG)) {
    fname=paste(prefix, REG$V4[r], suffix, sep="")
    DATA = read.table(fname, head=FALSE)
    cat (fname, "\n")		
    for (l in 1:nrow(DATA)) {
      A[DATA$V1[l]] = A[DATA$V1[l]] + DATA$V2[l]
      B[DATA$V1[l]] = B[DATA$V1[l]] + DATA$V3[l]
    }
  }
  bins=cut(C, breaks=BIN, labels=1:(nBIN-1))
  merged = cbind(A, B, bins)
  aggrE = as.vector(by(merged[, 1], merged[, 3], sum))
  aggrT = as.vector(by(merged[, 2], merged[, 3], sum))
  ser = aggrE*100/aggrT
  return (ser);
}

#LOAD SER BETWEEN ALL VARIANTS [ALWAYS *.fqt.* files here!!!!]
prefix="../../Source_data/WGS/Beagle5.4/N147754/benchmark_ukb23352_c20_qc_v1.subset.N147754.fullchr.";
suffix=".fqt.frequency.switch.txt.gz"
BGLser0=freqSER0(prefix, suffix, 147754);

prefix="../../Source_data/WGS/Shapeit5/N147754/benchmark_ukb23352_c20_qc_v1.subset.N147754.fullchr.shapeit5.ligated.";
suffix=".fqt.frequency.switch.txt.gz"
SHPser0=freqSER0(prefix, suffix, 147754);


jpeg("../../Supplementary_Information/Supplementary_figure3.jpeg", 3500, 1700, quality = 100, res=300)

X=1:(nBIN-1)
ser_bin=c(0.1, 0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0)
par(mfrow=c(1,2))
plot(X, log10(BGLser0), type="n", pch=20, xlab="Minor Allele Count", ylab="Switch Error Rate (%)", col="black", 
     lwd=2, xaxt='n', yaxt='n', main=paste("Switch Error Rate\n[N=147,754]", sep=""))
put.fig.letter(label="a", cex=2)
#abline(h=log10(ser_bin), col="lightgrey", lty=2)
#abline(v=X, col="lightgrey", lty=2)
text(X, par("usr")[3], labels = lBIN[X], srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.75)
axis(2, at=log10(ser_bin), label=ser_bin, las=2)
points(X, log10(BGLser0), type="o", pch=20, col="black", lwd=2, xaxt='n', yaxt='n')
points(X, log10(SHPser0), type="o", pch=20, col=COLdiff[2], lwd=2)
legend("topright", fill=c(COLdiff[2], "black"), legend=c("SHAPEIT5", "Beagle5.4"), bg="white")

plot(X, -((SHPser0 - BGLser0) * 100.0 / BGLser0), type="n", pch=20, xlab="Minor Allele Count",
     ylab = "Switch Error Rate Reduction (%)", col="black", lwd=2, xaxt='n', ylim=c(0, 65), main=paste("SER reduction\n[N=147,754]", sep=""), yaxt='n')
put.fig.letter(label="b", cex=2)
#abline(h=seq(-70, 70, 5), col="lightgrey", lty=2)
#abline(v=X, col="lightgrey", lty=2)
text(X, par("usr")[3], labels = lBIN[X], srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.75) 
axis(2, at=seq(-10,70,10), label=seq(-10,70,10), las=2)
points(X, -((SHPser0 - BGLser0) * 100.0 / BGLser0), type="o", pch=20, col=COLdiff[2], lwd=2)
dev.off()

