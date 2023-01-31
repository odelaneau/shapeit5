library(RColorBrewer)

pdf("FigureN20000.pdf", 7.5, 5)

#########################################################################################
#					Helper 						#
#					funtions			 		#
#########################################################################################

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


#########################################################################################
#					R2 imputation vs MAC				#
#########################################################################################
COLpair = brewer.pal(12,"Paired")
COLdiff = brewer.pal(8,"Set1")
COLsp = brewer.pal(10,"Spectral")
COLblues = brewer.pal(11,"RdBu")
COLp = c("black", COLdiff[2])

BIN=c(0,1,5,10,20,50,100,200,500,1000,2000,5000,10000,20000)
lBIN=c("singleton","2-5","6-10","11-20","21-50","51-100","101-200","201-500","501-1k","1k-2k","2k-5k","5k-10k","10k-20k")
nBIN=length(BIN)
options(scipen=999)
	
plot(1,1,type="n", ylim=c(0.0,1),xlim=c(1,(nBIN-1)), ylab="", xlab="UK Biobank Reference Panel Minor Allele Count", xaxt='n', las=1, main="WGS - Imputation accuracy")
mtext("UKB Reference panel N=19,881 | Target N=1,000", side=3, line=0.1)
mtext(expression(paste("Aggregate ", italic("r"),""^"2"*"")), side=2, line=2.2, cex=0.8)
text(1:(nBIN-1), par("usr")[3], labels = lBIN, srt = 45, adj = c(1.1,1.1), xpd = TRUE) 
abline(h=seq(0, 1, 0.2), col="lightgrey", lty=3)
abline(v=1:(nBIN-1), col="lightgrey", lty=3)

r=1
file=paste("/home/srubinac/fuse/PhasingUKB/Phasing/PhasingWGS/step8_imputation/reference_panels_N20000/concordance_test/binning2/imputed_1k_rp_beagle5.4_chr20_rp_af_binning2.rsquare.grp.txt.gz", sep="")
D = read.table(file, head=FALSE, stringsAsFactors=FALSE)
points(1:(nBIN-1), D$V5, type="o", pch=20, col=COLp[r], lwd=2)
r2_wgs_b5 <- D$V5
r=r+1
file=paste("/home/srubinac/fuse/PhasingUKB/Phasing/PhasingWGS/step8_imputation/reference_panels_N20000/concordance_test/binning2/imputed_1k_rp_shapeit5_chr20_rp_af_binning2.rsquare.grp.txt.gz", sep="")
D = read.table(file, head=FALSE, stringsAsFactors=FALSE)
points(1:(nBIN-1), D$V5, type="o", pch=20, col=COLp[r], lwd=2)
r2_wgs_s5 <- D$V5

rect(9.45, 0.11, 16.5,0.725, col="white",border = NA) #for subplot

legend("topleft", legend=c("UKB - SHAPEIT5.0","UKB - Beagle5.4"), fill=c(COLdiff[2],"black"), title="Reference panel", bg="white")

#subplot C
par(fig=c(0.24,0.497,0.05,0.4),new=TRUE, mar=c(5,4.2,4.2,2))
plot(1:(nBIN-1), (r2_wgs_s5 - r2_wgs_b5), type="n", pch=20, xlab="", ylab = "", col="black", lwd=2, xaxt='n', ylim=c(0, 0.1), xlim=c(1, 9), yaxt='n', bg="white")
mtext(expression(paste(italic("r"),""^"2"*"", " difference")), side=2, line=2.7, cex=0.8)
#abline(h=seq(0, 100, 5), col="lightgrey", lty=3)
#abline(v=1:(nBIN-1), col="lightgrey", lty=3)
text(1:9, par("usr")[3], labels = lBIN[1:9], srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.7) 
axis(2, at=c(0, 0.02, 0.04, 0.06, 0.08), las=2, cex=0.6)
points(1:9, (r2_wgs_s5[1:9]-r2_wgs_b5[1:9]), type="b", pch=20, col=COLdiff[3], lwd=2)
abline(h=0, col="black", lty=2)
title("Imputation accuracy increase",line=-1)

dev.off()

