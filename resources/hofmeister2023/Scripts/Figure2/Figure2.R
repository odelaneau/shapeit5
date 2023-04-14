setwd('~/Dropbox/SHAPEIT5/Final_revision/Figures/Scripts/Figure2/')

library(RColorBrewer)

pdf("../../Main_figures/Figure2.pdf", 15, 10)
par(mfrow=c(2,2), bg ="white")

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
#					FIGURE 2A					#
#					SER in WGS vs MAC		 		#
#########################################################################################

load(file = "../../Source_data/Main_fig2/data.Rdata")


COLpair = brewer.pal(12,"Paired")
COLdiff = brewer.pal(8,"Set1")

REG=read.table("../../Source_data/Main_fig2/chr20.size4Mb.txt", head=FALSE)
BIN=c(0,1,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000, 147754)
lBIN=c("singleton","2-5","6-10","11-20","21-50","51-100","101-200","201-500","501-1k","1k-2k","2k-5k","5k-10k","10k-20k","20k-50k","50k-100k", "100k+")
nREG=length(REG)
nBIN=length(BIN)
options(scipen=999)

# LOAD SER BETWEEN ALL VARIANTS
# prefix="../../Source_data/Main_fig2/Beagle5.4/N147754/benchmark_ukb23352_c20_qc_v1.subset.N147754.";
# suffix=".fqf.frequency.switch.txt.gz"
# 
# prefix="../../Source_data/Main_fig2/Shapeit5/N147754/benchmark_ukb23352_c20_qc_v1.subset.N147754.";
# suffix=".shapeit5.default.fqf.frequency.switch.txt.gz"

ser_bin=c(0.1, 0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0)
plot(log10(BGLser0), type="n", pch=20, xlab="Minor Allele Count", ylab="", col="black", lwd=2, xaxt='n', yaxt='n', ylim=log10(c(0.2,45)), main="UKB WGS - Phasing accuracy")
put.fig.letter(label="a", cex=2,font=2)
mtext("N=147,754", side=3, line=0.1)
mtext("Switch Error Rate (%)", side=2, line=2.2, cex=0.8)
#abline(h=log10(ser_bin), col="lightgrey", lty=3)
#abline(v=1:(nBIN-1), col="lightgrey", lty=3)
text(1:(nBIN-1), par("usr")[3], labels = lBIN, srt = 45, adj = c(1.1,1.1), xpd = TRUE) 
axis(2, at=log10(ser_bin), label=ser_bin, las=2)
points(log10(BGLser0), type="b", pch=20, col="black", lwd=2, xaxt='n', yaxt='n')
points(log10(SHPser0), type="b", pch=20, col=COLdiff[2], lwd=2)
log10.axis(2, at=seq(-1, 2)) ##REMOVE

legend("bottomleft", legend=c("SHAPEIT5.0","Beagle5.4"), fill=c(COLdiff[2],"black"), title="Phasing method", bg="white")

rect(9.45, log10(4.9), 16.5,log10(60), col="white",border = NA) #for subplot

#########################################################################################
#					FIGURE 2B					#
#					SER in WES vs MAC				#
#########################################################################################

library(RColorBrewer)
COLpair = brewer.pal(12,"Paired")
COLdiff = brewer.pal(8,"Set1")

BIN=c(0,1,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000, 500000)
lBIN=c("singleton","2-5","6-10","11-20","21-50","51-100","101-200","201-500","501-1k","1k-2k","2k-5k","5k-10k","10k-20k","20k-50k","50k-100k", "100k+")
nREG=length(REG)
nBIN=length(BIN)
options(scipen=999)

#LOAD SER BETWEEN ALL VARIANTS
# prefix="../Source_data/WES/Beagle5.4/SER.chr";
# suffix=".fq.frequency.switch.txt.gz"
# 
# prefix="../Source_data/WES/Shapeit5/SER_v2.chr";
# suffix=".cut0.001.prob_0.5.fq.frequency.switch.txt.gz"

ser_bin=c(0.1, 0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0)
plot(log10(BGLser1), type="n", pch=20, xlab="Minor Allele Count", ylab="", col="black", lwd=2, xaxt='n', yaxt='n', ylim=log10(c(0.2,45)), main="UKB WES - Phasing accuracy")
put.fig.letter(label="b", cex=2,font=2)
mtext("N=447,470", side=3, line=0.1)
mtext("Switch Error Rate (%)", side=2, line=2.2, cex=0.8)
#abline(h=log10(ser_bin), col="lightgrey", lty=3)
#abline(v=1:(nBIN-1), col="lightgrey", lty=3)
text(1:(nBIN-1), par("usr")[3], labels = lBIN, srt = 45, adj = c(1.1,1.1), xpd = TRUE) 
axis(2, at=log10(ser_bin), label=ser_bin, las=2)
points(log10(BGLser1), type="b", pch=20, col="black", lwd=2, xaxt='n', yaxt='n')
points(log10(SHPser1), type="b", pch=20, col=COLdiff[2], lwd=2)
log10.axis(2, at=seq(-1, 2)) ##REMOVE

legend("bottomleft", legend=c("SHAPEIT5.0","Beagle5.4"), fill=c(COLdiff[2],"black"), title="Phasing method", bg="white")

rect(9.45, log10(4.9), 16.5,log10(60), col="white",border = NA) #for subplot

#########################################################################################
#					FIGURE 2C					#
#					R2 imputation vs MAC				#
#########################################################################################
COLpair = brewer.pal(12,"Paired")
COLdiff = brewer.pal(8,"Set1")
COLsp = brewer.pal(10,"Spectral")
COLblues = brewer.pal(11,"RdBu")
COLp = c("black", COLdiff[2])

BIN=c(0,1,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,500000)
lBIN=c("singleton","2-5","6-10","11-20","21-50","51-100","101-200","201-500","501-1k","1k-2k","2k-5k","5k-10k","10k-20k","20k-50k","50k-100k", "100k+")
nBIN=length(BIN)
options(scipen=999)
	
plot(1,1,type="n", ylim=c(0.0,1),xlim=c(1,(nBIN-1)), ylab="", xlab="UK Biobank Reference Panel Minor Allele Count", xaxt='n', las=1, main="WGS - Imputation accuracy")
put.fig.letter(label="c", cex=2,font=2)
mtext("UKB Reference panel N=146,754 | Target N=1,000", side=3, line=0.1)
mtext(expression(paste("Aggregate ", italic("r"),""^"2"*"")), side=2, line=2.2, cex=0.8)
text(1:(nBIN-1), par("usr")[3], labels = lBIN, srt = 45, adj = c(1.1,1.1), xpd = TRUE) 
#abline(h=seq(0, 1, 0.2), col="lightgrey", lty=3)
#abline(v=1:(nBIN-1), col="lightgrey", lty=3)

file=paste("../../Source_data/Main_fig2/imputed_1k_rp_hrc_chr20_rp_af_binning2.rsquare.grp.txt.gz", sep="")
D = read.table(file, head=FALSE, stringsAsFactors=FALSE)
points(1:(nBIN-1), D$V5, type="o", pch=20, col="gray", lwd=2)

r=1
file=paste("../../Source_data/Main_fig2/imputed_1k_rp_beagle5.4_chr20_rp_af_binning2.rsquare.grp.txt.gz", sep="")
D = read.table(file, head=FALSE, stringsAsFactors=FALSE)
points(1:(nBIN-1), D$V5, type="o", pch=20, col=COLp[r], lwd=2)
r2_wgs_b5 <- D$V5
r=r+1
file=paste("../../Source_data/Main_fig2/imputed_1k_rp_shapeit5_chr20_rp_af_binning2.rsquare.grp.txt.gz", sep="")
D = read.table(file, head=FALSE, stringsAsFactors=FALSE)
points(1:(nBIN-1), D$V5, type="o", pch=20, col=COLp[r], lwd=2)
r2_wgs_s5 <- D$V5

rect(9.45, 0.11, 16.5,0.725, col="white",border = NA) #for subplot

legend("topleft", legend=c("UKB - SHAPEIT5.0","UKB - Beagle5.4", "HRC - N=27,165"), fill=c(COLdiff[2],"black","gray"), title="Reference panel", bg="white")
#legend(0.2, 0.8, legend=c("HRC - N=27,165"), fill=c("gray"), bg="white")

#########################################################################################
#					FIGURE 2D					#
#					R2 imputation vs MAC				#
#########################################################################################


BIN=c(0,1,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,147754)
lBIN=c("singleton","2-5","6-10","11-20","21-50","51-100","101-200","201-500","501-1k","1k-2k","2k-5k","5k-10k","10k-20k","20k-50k","50k-100k", "100k+")
nBIN=length(BIN)
options(scipen=999)

plot(1,1,type="n", ylim=c(0.0,1), ,xlim=c(1,(nBIN-1)), ylab="", xlab="UK Biobank Reference Panel Minor Allele Count", xaxt='n', las=1, main="WES - Imputation accuracy")
put.fig.letter(label="d", cex=2,font=2)
mtext("UKB Reference panel N=446,470 | Target N=1,000", side=3, line=0.1)
mtext(expression(paste("Aggregate ", italic("r"),""^"2"*"")), side=2, line=2.2, cex=0.8)
text(1:(nBIN-1), par("usr")[3], labels = lBIN, srt = 45, adj = c(1.1,1.1), xpd = TRUE) 
#abline(h=seq(0, 1, 0.2), col="lightgrey", lty=3)
#abline(v=1:(nBIN-1), col="lightgrey", lty=3)

r=1
file=paste("../../Source_data/Main_fig2/imputed_1k_rp_beagle5.4_allchrs_rp_af_binning2.rsquare.grp.txt.gz", sep="")
D = read.table(file, head=FALSE, stringsAsFactors=FALSE)
points(1:(nBIN-1), D$V5, type="o", pch=20, col=COLp[r], lwd=2)
r2_wes_b5 <- D$V5
r=r+1
file=paste("../../Source_data/Main_fig2/imputed_1k_rp_shapeit5_allchrs_rp_af_binning2.rsquare.grp.txt.gz", sep="")
D = read.table(file, head=FALSE, stringsAsFactors=FALSE)
points(1:(nBIN-1), D$V5, type="o", pch=20, col=COLp[r], lwd=2)
r2_wes_s5 <- D$V5

legend("topleft", legend=c("UKB - SHAPEIT5.0","UKB - Beagle5.4"), fill=c(COLdiff[2],"black"), title="Reference panel", bg="white")

rect(9.45, 0.11, 16.5,0.725, col="white",border = NA) #for subplot


lBIN=c("sing.","2-5","6-10","11-20","21-50","51-100","101-200","201-500","501-1k","1k-2k","2k-5k","5k-10k","10k-20k","20k-50k","50k-100k", "100k+")
#subplot A
par(fig=c(0.24,0.497,0.70,1),new=TRUE, mar=c(5,4.2,4.2,2))
plot(1:(nBIN-1), (SHPser0 - BGLser0) * 100.0 / BGLser1, type="n", pch=20, xlab="", ylab = "", col="black", lwd=2, xaxt='n', ylim=c(0, 53), xlim=c(1, 9), main="", yaxt='n')
mtext("SER difference (%)", side=2, line=2.2, cex=0.7)
#abline(h=seq(-50, 50, 5), col="lightgrey", lty=3)
#abline(v=1:(nBIN-1), col="lightgrey", lty=3)
text(1:9, par("usr")[3], labels = lBIN[1:9], srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.7) 
axis(2, at=seq(0,50,10), label=seq(0,50,10), las=2)
points(1:9, -(SHPser0[1:9] - BGLser0[1:9]) * 100.0 / BGLser0[1:9], type="b", pch=20, col=COLdiff[3], lwd=2)
abline(h=0, col="black", lty=2)
title("Phasing accuracy increase",line=-8)

#subplot B
par(fig=c(0.74,.997,0.70,1),new=TRUE, mar=c(5,4.2,4.2,2))
plot(1:(nBIN-1), (SHPser1 - BGLser1) * 100.0 / BGLser1, type="n", pch=20, xlab="", ylab = "", col="black", lwd=2, xaxt='n', ylim=c(0, 53), xlim=c(1, 9), main="", yaxt='n')
mtext("SER difference (%)", side=2, line=2.2, cex=0.7)
#abline(h=seq(-50, 50, 5), col="lightgrey", lty=3)
#abline(v=1:(nBIN-1), col="lightgrey", lty=3)
text(1:9, par("usr")[3], labels = lBIN[1:9], srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.7) 
axis(2, at=seq(0,50,10), label=seq(0,50,10), las=2)
points(1:9, -(SHPser1[1:9] - BGLser1[1:9]) * 100.0 / BGLser1[1:9], type="b", pch=20, col=COLdiff[4], lwd=2)
abline(h=0, col="black", lty=2)
title("Phasing accuracy increase",line=-8)

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

#subplot D
par(fig=c(0.74,0.997,0.05,0.4),new=TRUE, mar=c(5,4.2,4.2,2))
plot(1:(nBIN-1), (r2_wes_s5 - r2_wes_b5), type="n", pch=20, xlab="", ylab = "", col="black", lwd=2, xaxt='n', ylim=c(0, 0.1), xlim=c(1, 9), yaxt='n', bg="white")
mtext(expression(paste(italic("r"),""^"2"*"", " difference")), side=2, line=2.7, cex=0.8)
#abline(h=seq(0, 100, 5), col="lightgrey", lty=3)
#abline(v=1:(nBIN-1), col="lightgrey", lty=3)
text(1:9, par("usr")[3], labels = lBIN[1:9], srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.7) 
axis(2, at=c(0, 0.02, 0.04, 0.06, 0.08, 0.1), las=2, cex=0.6)
points(1:9, (r2_wes_s5[1:9]-r2_wes_b5[1:9]), type="b", pch=20, col=COLdiff[4], lwd=2)
abline(h=0, col="black", lty=2)
title("Imputation accuracy increase",line=-1)


dev.off()

