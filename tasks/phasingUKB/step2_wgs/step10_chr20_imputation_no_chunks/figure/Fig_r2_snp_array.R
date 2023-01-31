#DATA LOADING
library(RColorBrewer)
COLdiff = brewer.pal(8,"Set1")

CONFMall = read.csv("methods_coding.csv", stringsAsFactors=FALSE)

mets= c("BEAGLE_v5.4","SHAPEIT_v5.0.0")
refs = c("149,119 WGS samples")

list_legend=""
list_colours=""
j=1
	
	pdf(paste("Fig_r2_snp_array.pdf", sep=""), 6, 6, bg="white")
	is_new=FALSE
	par(fig=c(0,1,0,1), new=is_new, mar=c(5,5,4,2)+0.1)
	plot(1,1,type="n", xlim=c(0.0005,50), ylim=c(0.0,1), ylab=expression(paste("Aggregate ", italic("r"),""^"2"*"")), xlab="British Population Minor Allele Frequency (%)", log="x", xaxt='n', las=1)
	mtext(paste0("Target: 1000 samples | Reference:", refs[j]), side=3, line=0)
	axis(1, at=c(0.0001, 0.001, 0.01, 0.1, 1, 10, 50), label=c("0.0001", "0.001","0.01", "0.1", "1", "10", "50"), main=paste0("Imputation performance"))
	abline(v=c(0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50), col="lightgrey", lty=3)
	abline(h=seq(0,1,0.05), col="lightgrey", lty=3)
	x0=0


		
	listN<-c()
	listC<-c()
	#for (i in seq(1, 2, by=1)) 
	#{
		SUBmetP <- paste(mets,"_1k",sep="")
		CONFmetP = CONFMall[which(CONFMall$Prefix %in% SUBmetP),]
		for (r in 1:nrow(CONFmetP)) 
		{
		  file=paste("data/imputed_1k_rp_",CONFmetP$Met1[r],"_chr20.rsquare.grp.txt.gz", sep="")
		  print(file)
		  if(!file.exists(file))
		  {
			next
		  }
		  D = read.table(file, head=FALSE, stringsAsFactors=FALSE)
		  points(D$V3*100, D$V5, type="o", pch=20, col=COLdiff[r], lwd=2)

		}
		listN<-c(listN, paste0(CONFmetP$Met, "-phased reference panel"))
		listC<-c(listC, COLdiff[1:3])
	#}
	print(listN)
	print(listC)
	legend("bottomright", legend=listN, fill=listC, bg="white",cex = 0.8)
	
	dev.off()
#}
