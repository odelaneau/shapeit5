#########################################################################################
#		LOAD DATA FOR FIGURE 2A							#
#########################################################################################

REG=read.table("../../Data/WGS/chr20.size4Mb.txt", head=FALSE)
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

#LOAD SER BETWEEN ALL VARIANTS [ALWAYS *.fqf.* files here!!!!]
prefix="../../Data/WGS/Beagle5.4/N147754/benchmark_ukb23352_c20_qc_v1.subset.N147754.fullchr.";
suffix=".fqf.frequency.switch.txt.gz"
BGLser0=freqSER0(prefix, suffix, 147754);

prefix="../../Data/WGS/Shapeit5/N147754/benchmark_ukb23352_c20_qc_v1.subset.N147754.fullchr.shapeit5.ligated.";
suffix=".fqf.frequency.switch.txt.gz"
SHPser0=freqSER0(prefix, suffix, 147754);


#########################################################################################
#		LOAD DATA FOR FIGURE 2B							#
#########################################################################################

freqSER1 <- function (prefix, suffix, n) {
	A = rep(0, n)
	B = rep(0, n)
	C = 1:n
	for (r in 1:22) {
		fname=paste(prefix, r, suffix, sep="")
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

#LOAD SER BETWEEN ALL VARIANTS
prefix="../../Data/WES/Beagle5.4/Trios_only/benchmark_beagle_Trios.chr";
suffix=".fq.frequency.switch.txt.gz"
BGLser1=freqSER1(prefix, suffix, 500000);

prefix="../../Data/WES/Shapeit5/Trios_only/benchmark_Trios.chr";
suffix=".cut0.001.prob_0.5.fq.frequency.switch.txt.gz"
SHPser1=freqSER1(prefix, suffix, 500000);


save(list = c("SHPser0", "SHPser1","BGLser0", "BGLser1"), file = "data.Rdata")

