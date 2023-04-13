setwd('~/Dropbox/SHAPEIT5/Final_revision/Figures/Scripts/extended_figure6/')
# 02-Sep-2022 Diogo Ribeiro @ UNIL
# Basic stats on LoF variants

library(data.table)
library(ggplot2)

data = fread("../../Source_data/CH/lof.results.vars.txt.gz", header = T)
indsData = fread("../../Source_data/CH/lof.results.inds.txt.gz", header = T)
totalInds = 374826

freq = data.table(table(data$gene))

jpeg("../../Extended_figures/Extended_figure6.jpeg", 3000, 1000, quality = 100, res=300)

## PANEL a
g1 = ggplot(freq,aes(N)) + 
  geom_histogram(bins = 50, color = "white") +
  annotate("text",label = paste("Median = ", median(freq$N), "\nMean = ", round(mean(freq$N),1)), x = 100, y = 1200, fontface = "bold", size = 3) +
  theme_minimal() +
  ggtitle("a. # LoF variants per gene") +
  xlab("# LoF variants") +
  ylab("# Genes") +
  # labs(tag = "a") + 
  xlim(c(0,150)) + 
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5, size = 12, face = "bold"), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1  )

## PANEL b
## plot AF (categories/bins) per LoF variant
data$ac = data$af * totalInds * 2
data$ac_cat = "-1"
data[ac == 1]$ac_cat = "singl."
data[ac == 2]$ac_cat = "doubl."
data[ac >= 3][ac <= 5]$ac_cat = "3-5"
data[ac > 5][ac <= 10]$ac_cat = "6-10"
data[ac > 10][ac <= 20]$ac_cat = "10-20"
data[ac > 20][ac <= 50]$ac_cat = "20-50"
data[ac > 50][ac <= 100]$ac_cat = "50-100"
data[ac > 100]$ac_cat = ">100"

acData = data.table(table(data$ac_cat))
g2 = ggplot(acData,aes(x = V1, y = N)) + 
  geom_bar(stat = "identity", fill = "grey", color = "black", size = 0.5) +
  scale_x_discrete(limits = c("singl.","doubl.","3-5","6-10","10-20","20-50","50-100",">100")) + 
  theme_minimal() +
  ggtitle("b. # LoF per allele count") +
  xlab("Allele count") +
  ylab("# LoF variants") +
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5, size = 12, face = "bold"), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1, axis.text.x = element_text(angle = 20, vjust=0.6)  )

## PANEL c
## Per individual
freqInds = data.table(table(indsData$individual))
rest = data.table(V1 = "missing", N = rep(0,totalInds - nrow(freqInds)))
freqInds = rbind(freqInds,rest)
freqFreq = data.table(table(freqInds$N))
g3 = ggplot(freqFreq,aes(x = as.numeric(V1),y = N)) + 
  geom_bar( stat = "identity") +
  annotate("text",label = paste("Median = ", median(freqInds$N), "\nMean = ", round(mean(freqInds$N),1)), x = 15, y = 40000, fontface = "bold", size = 3) +
  xlim(c(0,20)) +
  ggtitle("c. # LoF per individual") +
  xlab("# LoF variants") +
  ylab("# individuals") +
  theme_minimal() +
  theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5, size = 12, face = "bold"), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1  )

g1 + g2 + g3
dev.off()
