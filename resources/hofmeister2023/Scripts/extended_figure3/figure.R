###### COLORS ######
setwd('~/Dropbox/SHAPEIT5/Final_revision/Figures/Scripts/extended_figure3/')
library(RColorBrewer)
COLpair = brewer.pal(12,"Paired")
COLdiff = brewer.pal(8,"Set1")
COLlabel= c(rgb(0.3,0.9,0.4,0.6), rgb(0.3,0.5,0.4,0.6), rgb(0.3,0.1,0.4,0.6))

Ntrio = read.table("../../Source_data/Singleton/dataWGS_trios.txt", head=FALSE)
Nduo = read.table("../../Source_data/Singleton/dataWGS_duos.txt", head=FALSE)
SER = read.table("../../Source_data/Singleton/dataWGS_ser.txt", head=FALSE)


#pdf("figure.pdf", 8,8)
jpeg("../../Extended_figures/Extended_figure3.jpeg", 3000, 3000, quality = 100, res=300)
par(mfrow=c(2,2))

Y = Nduo$V1*100/sum(Nduo$V1)
myBar = barplot(Y,space = 0.5, border=F , xlim=c(0,5), ylim=c(0,60), col=COLlabel, ylab="Percentage (%)", names.arg=c("p=0/0", "p=0/1", "p=1/1"), xlab="Parental genotype when offspring is 0/1", main="a. Singleton phasing using Duos\nWGS data - N=147,754")
text(myBar, Y+2 , paste("n=", Nduo$V1, sep="") ,cex=1) 

Y = Ntrio$V1*100/sum(Ntrio$V1)
myBar = barplot(Y,space = 0.5, border=F , xlim=c(0,5), ylim=c(0,60), col=COLlabel, ylab="Percentage (%)", names.arg=c("f=0/1\nm=0/0", "f=0/0\nm=0/1", "f=0/0\nm=0/0"), xlab="Parental genotypes when offspring is 0/1", main="b. Singleton phasing using Trios\nWGS data - N=147,754")
text(myBar, Y+2 , paste("n=", Ntrio$V1, sep="") ,cex=1) 

Y=as.vector(c(SER[1,2], SER[1,3]-SER[1,2]))
Y = Y*100/sum(Y)
myBar = barplot(Y,space = 0.5, border=F , xlim=c(0,3.5), ylim=c(0,80), col=COLlabel, ylab="Switch error rate (%)", names.arg=c("Incorrect", "Correct"), xlab="Singleton phasing", main="c. Singleton phasing using SHAPEIT5\nValidation from Duos\nWGS data - N=147,754")
abline(h=50)
text(myBar, Y+2 , paste("n=", c(SER[1, 2], SER[1,3]-SER[1,2]), sep="") ,cex=1) 

legend("topleft", bty="n", legend=paste0("binomial p-value\n=",formatC(binom.test(SER[1, 2],(SER[1, 2]+(SER[1,3]-SER[1,2])), alternative="two.sided")$p.value, format = "e", digits = 2)))

Y=as.vector(c(SER[2,2], SER[2,3]-SER[2,2]))
Y = Y*100/sum(Y)
myBar = barplot(Y,space = 0.5, border=F , xlim=c(0,3.5), ylim=c(0,80), col=COLlabel, ylab="Switch error rate (%)", names.arg=c("Incorrect", "Correct"), xlab="Singleton phasing", main="d. Singleton phasing using SHAPEIT5\nValidation from Trios\nWGS data - N=147,754")
abline(h=50)
text(myBar, Y+2 , paste("n=", c(SER[2, 2], SER[2,3]-SER[2,2]), sep="") ,cex=1) 
legend("topleft", bty="n", legend=paste0("binomial p-value\n=",formatC(binom.test(SER[2, 2],(SER[2, 2]+(SER[2,3]-SER[2,2])), alternative="two.sided")$p.value, format = "e", digits = 2)))


dev.off()

