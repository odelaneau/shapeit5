#DATA LOADING
library(RColorBrewer)
COLdiff = brewer.pal(8,"Set1")

CONFMall = read.csv("methods_coding.csv", stringsAsFactors=FALSE)

mets= c("BEAGLE_v5.3","SHAPEIT_v5.0.0","SHAPEIT_v5.0.0_robin")
refs = c("147654 UKB WGS samples")

list_legend=""
list_colours=""
j=1
	pdf(paste("Fig_r2_snp_array_robin.pdf", sep=""), 6, 6, bg="white")
	plot(1,1,type="n", xlim=c(0.00008,50), ylim=c(0.0,1), ylab=expression(paste(italic("r"),""^"2"*" between imputed and 30x genotypes")), xlab="British Population Minor Allele Frequency (%)", log="x", xaxt='n', las=1)
	mtext(paste0("Target: 100 samples | Reference:", refs[j]), side=3, line=0)
	axis(1, at=c(0.0001, 0.001, 0.01, 0.1, 1, 10, 50), label=c("0.0001", "0.001","0.01", "0.1", "1", "10", "50"), main=paste0("Imputation performance"))
	abline(v=c(0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50), col="lightgrey", lty=3)
	abline(h=seq(0,1,0.05), col="lightgrey", lty=3)
	x0=0
	is_new=FALSE
	par(fig=c(0,1,0,1), new=is_new, mar=c(5,4.2,4,2)+0.1)
		
	listN<-c()
	listC<-c()
	#for (i in seq(1, 2, by=1)) 
	#{
		SUBmetP <- paste(mets,"_100",sep="")
		CONFmetP = CONFMall[which(CONFMall$Prefix %in% SUBmetP),]
		for (r in 1:nrow(CONFmetP)) 
		{
		  file=paste("data/imputed_100_rp_",CONFmetP$Met1[r],".rsquare.grp.txt.gz", sep="")
		  print(file)
		  if(!file.exists(file))
		  {
			next
		  }
		  D = read.table(file, head=FALSE, stringsAsFactors=FALSE)
		  points(D$V3*100, D$V5, type="o", pch=20, col=COLdiff[r], lwd=2)

		}
		listN<-c(listN, paste(CONFmetP$Met, " phased reference panel"))
		listC<-c(listC, COLdiff[1:3])
	#}
	print(listN)
	print(listC)
	legend("bottomright", legend=listN, fill=listC, bg="white",cex = 0.8)
	
	dev.off()
#}
