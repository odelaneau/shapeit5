setwd('~/Dropbox/SHAPEIT5/Final_revision/Figures/Scripts/extended_figure5')
library(ggplot2)
library(RColorBrewer)
COLpair = brewer.pal(12,"Paired")
COLdiff = brewer.pal(8,"Set1")
COLsp = brewer.pal(10,"Spectral")
COLblues = brewer.pal(11,"RdBu")
CB= c(brewer.pal(11,"RdYlBu")[1:5],brewer.pal(11,"RdYlBu")[8:11]) 

put.fig.letter <- function(label, location="topleft", x=NULL, y=NULL, 
                           offset=c(0, 0), ...) {
  if(length(label) > 1) {
    warning("length(label) > 1, using label[1]")
  }
  if(is.null(x) | is.null(y)) {
    coords <- switch(location,
                     topleft = c(0.015,0.98),
                     topcenter = c(0.5525,0.98),
                     topright = c(0.985, 0.98),
                     bottomleft = c(0.015, 0.02), 
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

log10.axis <- function(side, at, ...) {
    at.minor <- log10(outer(1:9, 10^(min(at):max(at))))
    lab <- sapply(at, function(i) as.expression(bquote(10^ .(i))))
    axis(side=side, at=at.minor, labels=NA, tcl=par("tcl")*0.5, ...)
    #axis(side=side, at=at, labels=lab, ...)
}


#PLOT1
f2<-read.table('../../Source_data/WGS/Summary_data/Summary_WGS_f2.txt', hea=T)
S<-read.table('../../Source_data/WGS/Summary_data/Summary_WGS_S.txt', hea=T)
f2_copy<-f2
jpeg("../../Extended_figures/Extended_figure5.jpeg", 3000, 3000, quality = 100, res=300)
par(mfrow=c(2,2))
f2<-f2[f2$prob!=1,]
f2<-f2[f2$bin<10,]
f2<-f2[f2$bin>1,]
probs<-c(0.5,0.6,0.7,0.75,0.8,0.85,0.9,0.95,0.99,1)
lBIN=c("2-5","6-10","11-20","21-50","51-100","101-200","201-500","501-1k","1k-2k","2k-5k","5k-10k","10k-20k","20k-50k","50k-100k", "100k-150k", "150k-250k", "250k-500k")
s=8
X = 1:(9+s)
X=1:(length(unique(f2$bin)))
ser_bin=c(1.0,2.0,5.0,10.0,15,20.0)
plot(rep(X, length(unique(f2$prob))), log10(f2$ser), type="n", pch=20, xlab="Minor Allele Count", ylab="", 
     col="black", lwd=2, xaxt='n', yaxt='n',
     main="WGS - accuracy using phasing probability filter\n(N=147,754)", xlim=c(1, 8 ), ylim=log10(c(1,15)))
put.fig.letter(label="a", cex=2)
mtext("Switch Error Rate (%)", side=2, line=2.2, cex=0.8)
#abline(h=log10(ser_bin), col="lightgrey", lty=3)
#abline(v=X, col="lightgrey", lty=3)
text(X, par("usr")[3], labels = lBIN[X], srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.75)
axis(2, at=log10(ser_bin), label=ser_bin, las=2)
log10.axis(2, at=seq(-1, 2)) ##REMOVE

C=0
for (PROB in probs[-length(probs)]){
  C=C+1
  points(X, log10(f2$ser[f2$prob==PROB]), type="o", pch=20, col=CB[C], lwd=2)
  
}
legend("topright", fill=CB[1:C], legend=c(probs[-length(probs)]), title="Filter", bg="white", cex=1)


#PLOT2

S<-read.table('../../Source_data/WGS/Summary_data/Summary_WGS_S.txt', hea=T)
S<-S[S$bin>1,]
X=1:(length(unique(S$bin)))
plot(X, S$prop[S$prob==0.5]*100, type="n", pch=20, xlab="Minor Allele Count", ylab="", col="black",
     lwd=2, xaxt='n',  yaxt='n', main=paste("WGS - missigness using phasing probability filter\n(N=147,754) "),  xlim=c(1, 8 ), ylim=c(0,25)) #, xlim=c(min(X), max(X)), ylim=c(min(S$prop*100), max(S$prob*100)))
put.fig.letter(label="b", cex=2)
mtext("Heterozygous sites removed after filtering (%)", side=2, line=2.2, cex=0.8)
#abline(h=seq(0,100,5), col="lightgrey", lty=3)
#abline(v=X, col="lightgrey", lty=3)
text(X, par("usr")[3], labels = lBIN[X], srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.75)
axis(2, at=seq(0,100,5), label=seq(0,100,5), las=2)
C=0
for (BIN in unique(S$prob)){
  C=C+1
  points(X, 100-S$prop[S$prob==BIN]*100, type="o", pch=20, col=CB[C], lwd=2)
}
legend("topright", fill=CB[1:C], legend=c(probs[-length(probs)]), title="Filter", bg="white", cex=1)


#PLOT3
f2<-read.table('../../Source_data/WES/Summary_data/Summary_WES_f2.txt', hea=T)
S<-read.table('../../Source_data/WES/Summary_data/Summary_WES_S.txt', hea=T)
f2_copy<-f2
f2<-f2[f2$prob!=1,]
f2<-f2[f2$bin<10,]
f2<-f2[f2$bin>1,]
probs<-c(0.5,0.6,0.7,0.75,0.8,0.85,0.9,0.95,0.99,1)
lBIN=c("2-5","6-10","11-20","21-50","51-100","101-200","201-500","501-1k","1k-2k","2k-5k","5k-10k","10k-20k","20k-50k","50k-100k", "100k-150k", "150k-250k", "250k-500k")
s=8
X = 1:(9+s)
X=1:(length(unique(f2$bin)))
ser_bin=c(0.5,1.0,2.0,5.0,10.0,20.0)
plot(rep(X, length(unique(f2$prob))), log10(f2$ser), type="n", pch=20, xlab="Minor Allele Count", ylab="", 
     col="black", lwd=2, xaxt='n', yaxt='n',
     main="WES - accuracy using phasing probability filter\n(N=447,470)", xlim=c(1, 8 ), ylim=log10(c(0.5,20.5)))
put.fig.letter(label="c", cex=2)
mtext("Switch Error Rate (%)", side=2, line=2.2, cex=0.8)
#abline(h=log10(ser_bin), col="lightgrey", lty=3)
#abline(v=X, col="lightgrey", lty=3)
text(X, par("usr")[3], labels = lBIN[X], srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.75)
axis(2, at=log10(ser_bin), label=ser_bin, las=2)
log10.axis(2, at=seq(-1, 2)) ##REMOVE

C=0
for (PROB in probs[-length(probs)]){
  C=C+1
  points(X, log10(f2$ser[f2$prob==PROB]), type="o", pch=20, col=CB[C], lwd=2)
  
}
legend("topright", fill=CB[1:C], legend=c(probs[-length(probs)]), title="Filter", bg="white", cex=1)


#PLOT4

S<-read.table('../../Data/WES/Summary_data/Summary_WES_S.txt', hea=T)
S<-S[S$bin>1,]
X=1:(length(unique(S$bin)))
plot(X, S$prop[S$prob==0.5]*100, type="n", pch=20, xlab="Minor Allele Count", ylab="", col="black",
     lwd=2, xaxt='n',  yaxt='n', main=paste("WES - missigness using phasing probability filter\n(N=447,470) "),  xlim=c(1, 8 ), ylim=c(0,58)) #, xlim=c(min(X), max(X)), ylim=c(min(S$prop*100), max(S$prob*100)))
put.fig.letter(label="d", cex=2)
mtext("Heterozygous sites removed after filtering (%)", side=2, line=2.2, cex=0.8)
#abline(h=seq(0,100,5), col="lightgrey", lty=3)
#abline(v=X, col="lightgrey", lty=3)
text(X, par("usr")[3], labels = lBIN[X], srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.75)
axis(2, at=seq(0,100,5), label=seq(0,100,5), las=2)
C=0
for (BIN in unique(S$prob)){
  C=C+1
  points(X, 100-S$prop[S$prob==BIN]*100, type="o", pch=20, col=CB[C], lwd=2)
}
legend("topright", fill=CB[1:C], legend=c(probs[-length(probs)]), title="Filter", bg="white", cex=1)


dev.off()





