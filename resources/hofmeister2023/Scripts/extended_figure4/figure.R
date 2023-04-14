setwd('~/Dropbox/SHAPEIT5/Final_revision/Figures/Scripts/extended_figure4/')
library(ggplot2)
library(RColorBrewer)
library(data.table)
COLpair = brewer.pal(12,"Paired")
COLdiff = brewer.pal(8,"Set1")
COLsp = brewer.pal(10,"Spectral")
COLblues = brewer.pal(11,"RdBu")

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

jpeg("../../Extended_figures/Extended_figure4.jpeg", 3000, 3000, quality = 100, res=300)
par(mfrow=c(2,1), bg ="white")

COLpair = brewer.pal(12,"Paired")
COLdiff = brewer.pal(8,"Set1")
COLsp = brewer.pal(10,"Spectral")
COLblues = brewer.pal(11,"RdBu")
COLp = c("black", COLdiff[2])

BIN=c(0,1,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,500000)
lBIN=c("singleton","2-5","6-10","11-20","21-50","51-100","101-200","201-500","501-1k","1k-2k","2k-5k","5k-10k","10k-20k","20k-50k","50k-100k", "100k+")
nBIN=length(BIN)
options(scipen=999)

file=paste("../../Source_data/Extended_fig4/wgs_nrd.txt", sep="")
D = read.table(file, head=TRUE, stringsAsFactors=FALSE)
transpose(D)
D2 = df = subset(D, select = -c(bin) )	
DT <- transpose(D2)
MT <- as.matrix(data.frame(DT))
colnames(MT) <- c(lBIN)

x <- barplot(MT,beside = TRUE, col=c(COLdiff[2],"black","gray"), ylim=c(0,100),las=2,main="WGS - Imputation accuracy",xlab="UK Biobank Reference Panel Minor Allele Count",xaxt="n")
put.fig.letter(label="a", cex=2)
text(x=x[2,]-.25, par("usr")[3], lBIN, xpd=TRUE, srt=45, adj = c(1.1,1.1))
mtext("UKB Reference panel N=146,754 | Target N=1,000", side=3, line=0.1)
mtext("Non-reference discordance rate %", side=2, line=2.2)

legend("topright", legend=c("UKB - SHAPEIT5.0","UKB - Beagle5.4", "HRC - N=27,165"), fill=c(COLdiff[2],"black","gray"), title="Reference panel", bg="white")

#########################################################################################
#					FIGURE 2D					#
#					R2 imputation vs MAC				#
#########################################################################################


BIN=c(0,1,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,147754)
lBIN=c("singleton","2-5","6-10","11-20","21-50","51-100","101-200","201-500","501-1k","1k-2k","2k-5k","5k-10k","10k-20k","20k-50k","50k-100k", "100k+")
nBIN=length(BIN)
options(scipen=999)

file=paste("../../Source_data/Extended_fig4/wes_nrd.txt", sep="")
D = read.table(file, head=TRUE, stringsAsFactors=FALSE)
transpose(D)
D2 = df = subset(D, select = -c(bin) )	
DT <- transpose(D2)
MT <- as.matrix(data.frame(DT))
colnames(MT) <- c(lBIN)

x <- barplot(MT,beside = TRUE, col=c(COLdiff[2],"black"), ylim=c(0,100),las=2,main="WES - Imputation accuracy",xlab="UK Biobank Reference Panel Minor Allele Count",xaxt="n")
put.fig.letter(label="b", cex=2)
text(x=x[2,]-.25, par("usr")[3], lBIN, xpd=TRUE, srt=45, adj = c(1.1,1.1))
mtext("UKB Reference panel N=446,470 | Target N=1,000", side=3, line=0.1)
mtext("Non-reference discordance rate %", side=2, line=2.2)

legend("topright", legend=c("UKB - SHAPEIT5.0","UKB - Beagle5.4"), fill=c(COLdiff[2],"black"), title="Reference panel", bg="white")

dev.off()





