###### COLORS ######

library(RColorBrewer)
COLpair = brewer.pal(12,"Paired")
COLdiff = brewer.pal(8,"Set1")

###### VECTORS ######
REG=read.table("/home/olivier/Dropbox/Repository/shapeit5/tasks/phasingUKB/step2_wgs/step2_splitchunks/chr20.size4Mb.txt", head=FALSE)
SIZE=c(2000,5000,10000,20000,50000,100000,147754)
BIN=c(0,1,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000, 147754)
lBIN=c("singleton","2-5","6-10","11-20","21-50","51-100","101-200","201-500","501-1k","1k-2k","2k-5k","5k-10k","10k-20k","20k-50k","50k-100k", "100k+")
TYPE=c("fqa", "fqc", "fqr", "fqf")
nSIZE=length(SIZE)
nTYPE=length(TYPE)
nREG=length(REG)
nBIN=length(BIN)
options(scipen=999)

#########################################################################################
#					FIGURE2						#
#			Plot switch error rates	vs MAC	 				#						
#########################################################################################

freqSER <- function (prefix, midfix, suffix, s) {
	n = as.numeric(BIN[10+s])
	A = rep(0, n)
	B = rep(0, n)
	C = 1:n
	for (r in 1:nrow(REG)) {
		tmp=paste(prefix, midfix, REG$V4[r], suffix, sep=".");
		DATA = read.table(tmp, head=FALSE)
		cat (tmp, " ", nrow(DATA), " ", max(DATA$V1), "\n")		
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

freqSER2 <- function (fname, s) {
	n = as.numeric(BIN[10+s])
	A = rep(0, n)
	B = rep(0, n)
	C = 1:n
	DATA = read.table(fname, head=FALSE)
	cat (fname, " ", nrow(DATA), " ", max(DATA$V1), "\n")
	for (l in 1:nrow(DATA)) {
		A[DATA$V1[l]] = A[DATA$V1[l]] + DATA$V2[l]
		B[DATA$V1[l]] = B[DATA$V1[l]] + DATA$V3[l]
	}
	bins=cut(C, breaks=BIN, labels=1:(nBIN-1))
	merged = cbind(A, B, bins)
	aggrE = as.vector(by(merged[, 1], merged[, 3], sum))
	aggrT = as.vector(by(merged[, 2], merged[, 3], sum))
	ser = aggrE*100/aggrT
	return (ser);
}

#LOAD SER BETWEEN ALL VARIANTS
BQ0=rep(0, nSIZE)
BQ1=rep(0, nSIZE)
BQ2=rep(0, nSIZE)

#for (s in 1:nSIZE) {
for (s in nSIZE:nSIZE) {
	str1="~/Fuse/PhasingUKB/Phasing/PhasingWGS/step3_runbeagle/N";
	str2="/chunks/benchmark_ukb23352_c20_qc_v1.subset.N";
	prefix=paste(str1, SIZE[s], str2, SIZE[s], sep="")
	midfix="fullchr"
	suffix=paste(TYPE[4], ".frequency.switch.txt.gz", sep="")
	BGLser = freqSER(prefix, midfix, suffix, s);
			
	str1="~/Fuse/PhasingUKB/Phasing/PhasingWGS/step5_runshapeit_phase2/N";
	str2="/benchmark_ukb23352_c20_qc_v1.subset.N";
	prefix=paste(str1, SIZE[s], str2, SIZE[s], sep="")
	midfix="fullchr.shapeit5.ligated"
	suffix=paste(TYPE[4], ".frequency.switch.txt.gz", sep="")
	SHPser = freqSER(prefix, midfix, suffix, s);	
	
	str="~/Fuse/PhasingUKB/Phasing/PhasingWGS/step5_runshapeit_phase2/benchmark_ukb23352_c20_qc_v1.subset.N147754.fullchr.shapeit5.ligated.fqf.frequency.switch.txt.gz";
	SHPser2 = freqSER2(str, s);	
}



par(mfrow=c(1,2))
X=1:(nBIN-1)
ser_bin=c(0.1, 0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0)
plot(X, log10(BGLser), type="n", pch=20, xlab="Minor Allele Count", ylab="Switch Error Rate (%)", col="black", lwd=2, xaxt='n', yaxt='n', main=paste("A. Switch Error Rate\n[N=", SIZE[s], "]", sep=""))
abline(h=log10(ser_bin), col="lightgrey", lty=2)
abline(v=X, col="lightgrey", lty=2)
text(X, par("usr")[3], labels = lBIN[X], srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.75)
axis(2, at=log10(ser_bin), label=ser_bin, las=2)
points(X, log10(BGLser), type="o", pch=20, col="black", lwd=2, xaxt='n', yaxt='n')
points(X, log10(SHPser), type="o", pch=20, col=COLdiff[2], lwd=2)
points(X, log10(SHPser2), type="o", pch=20, col=COLdiff[3], lwd=2)
legend("topright", fill=c(COLdiff[3], COLdiff[2], "black"), legend=c("SHAPEIT5-fullchr", "SHAPEIT5-chunks", "Beagle5.4"), bg="white")

plot(X, (SHPser - BGLser) * 100.0 / BGLser, type="n", pch=20, xlab="Minor Allele Count", ylab = "Switch Error Rate Reduction (%)", col="black", lwd=2, xaxt='n', ylim=c(-50, 10), main=paste("B. SER reduction\n[N=", SIZE[s],"]", sep=""), yaxt='n')
abline(h=seq(-50, 50, 5), col="lightgrey", lty=2)
abline(v=X, col="lightgrey", lty=2)
text(X, par("usr")[3], labels = lBIN[X], srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.75) 
axis(2, at=seq(-60,10,10), label=seq(-60,10,10), las=2)
points(X, (SHPser - BGLser) * 100.0 / BGLser, type="o", pch=20, col=COLdiff[2], lwd=2)
points(X, (SHPser2 - BGLser) * 100.0 / BGLser, type="o", pch=20, col=COLdiff[3], lwd=2)
abline(h=0, col="black")




#########################################################################################
#					FIGURE2						#
#			Plot switch error rates at commons 				#						
#########################################################################################

for (s in nSIZE:nSIZE) {
	ERR0 = 0; TOT0 = 0;
	ERR1 = 0; TOT1 = 0;
	for (r in 1:nrow(REG)) {
		str=paste("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step3_runbeagle/N", SIZE[s], "/chunks/benchmark_ukb23352_c20_qc_v1.subset.N", SIZE[s], ".fullchr.", REG$V4[r], ".fqf.sample.switch.txt.gz", sep="")
		DATA = read.table(str, head=FALSE)
		cat (str, " ", nrow(DATA), "\n")
		DATA = DATA[which(!is.nan(DATA$V4)), ]
		ERR0 = ERR0 + sum(DATA$V2)
		TOT0 = TOT0 + sum(DATA$V3)
	}
	
	for (r in 1:nrow(REG)) {
		str=paste("~/Fuse/PhasingUKB/Phasing/PhasingWGS/step5_runshapeit_phase2/N", SIZE[s], "/benchmark_ukb23352_c20_qc_v1.subset.N", SIZE[s], ".fullchr.shapeit5.ligated.", REG$V4[r], ".fqf.sample.switch.txt.gz", sep="")
		DATA = read.table(str, head=FALSE)
		cat (str, " ", nrow(DATA), "\n")
		DATA = DATA[which(!is.nan(DATA$V4)), ]
		ERR1 = ERR1 + sum(DATA$V2)
		TOT1 = TOT1 + sum(DATA$V3)
	}
}



DATA=read.table("/home/olivier/Fuse/PhasingUKB/Phasing/PhasingWGS/step5_runshapeit_phase2/benchmark_ukb23352_c20_qc_v1.subset.N147754.fullchr.shapeit5.ligated.fqf.sample.switch.txt.gz", head=FALSE)

