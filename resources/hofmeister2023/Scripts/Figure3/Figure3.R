setwd('~/Dropbox/SHAPEIT5/Final_revision/Figures/Scripts/Figure3/')
#29-Sep-2022 Diogo Ribeiro @ UNIL
# CH results summary

library(data.table)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggsignif)

pdf("../../Main_figures/Figure3.pdf", 11, 8)

data1 = fread("../../Source_data/CH/lof.nosingle.results.txt")
data2 = fread("../../Source_data/CH/lof.nosingle.cutoff0.99.results.txt")
data3 = fread("../../Source_data/CH/lof.nosingle.rand.results.txt")

data <- fread("../../Source_data/CH/lof.enrichment.main.txt", stringsAsFactors = FALSE, header = TRUE, sep="\t")

dataLoF = fread("../../Source_data/CH/lof.nosingle.results.txt")
dataSyn = fread("../../Source_data/CH/synonymous.nosingle.results.txt")
dataMis = fread("../../Source_data/CH/missense.nosingle.results.txt")

labels = c("Full data","High confidence","Random phasing")
colors = c("#c6dbef","#2171b5","#fc9272")

#############
# Panel a
#############

n1CH = nrow(data1[numCH > 0])
n2CH = nrow(data2[numCH > 0])
n3CH = nrow(data3[numCH > 0])

dt = data.table(labels = labels, CH = c(n1CH,n2CH,n3CH), cols = colors)

g1 = ggplot(dt, aes(x = labels, y = CH, fill = cols)) + 
  geom_bar(stat = "identity", size = 0.7, color = "black", width = 0.7) +
  geom_text(aes(label = CH), size = 3, vjust = 1.5, fontface = "bold") +
  scale_fill_identity(labels = dt$labels, guide = guide_legend(reverse = TRUE)) + 
  scale_x_discrete(limits = dt$labels) +
  labs(tag = "a") +
  ylab("# compound het. genes") +
  xlab("Dataset") +
  theme_minimal() +
  theme(text = element_text(size = 12), plot.title = element_text(hjust = 0.5), 
        legend.position = "none",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill=NA, size = 1) )


#############
# Panel b
#############

options("scipen"=100, "digits"=2)

colnames(data) = c("ExternalList","Size1","ItemGroup","Size2","Overlap","OddsRatio","Pvalue")
MAX_X_AXIS = 6

data[ExternalList == "hart2014_essential"]$ExternalList = "Essential in culture\nHart 2014"
data[ExternalList == "hart2017_essential"]$ExternalList = "Essential CRISPR\nHart 2017"
data[ExternalList == "gnomad_essential"]$ExternalList = "Essential gnomAD\nKarczewski 2020"
data[ExternalList == "gnomad_nonessential"]$ExternalList = "Non-essential\ngnomAD\nKarczewski 2020"
data[ExternalList == "gnomad_homozygous_tolerant"]$ExternalList = "Homozygous\nLoF tolerant\nKarczewski 2020"
data[ExternalList == "mgi_essential"]$ExternalList = "Essential in mice\nGeorgi 2013"
data[ExternalList == "hart2014_nonessential"]$ExternalList = "Non-essential\nin culture\nHart 2014"
data[ExternalList == "ADaM_essential"]$ExternalList = "Essential ADaM\nVinceti 2021"

o = c("Homozygous\nLoF tolerant\nKarczewski 2020","Non-essential\ngnomAD\nKarczewski 2020","Non-essential\nin culture\nHart 2014",
      "Essential in culture\nHart 2014","Essential CRISPR\nHart 2017", "Essential in mice\nGeorgi 2013","Essential gnomAD\nKarczewski 2020","Essential ADaM\nVinceti 2021")
data$ExternalList = factor(as.character(data$ExternalList), levels = o)
data = data[order(-data$ExternalList)]

direct_fisher <- function(data, data1, data2) {
  dataset = data.table()
  for (extList in unique(data$ExternalList)){
    d1 = data[ExternalList == extList][ItemGroup == data1]
    d2 = data[ExternalList == extList][ItemGroup == data2]
    TP = d1$Overlap
    FP = d1$Size2 - d1$Overlap
    FN = d2$Overlap
    TN = d2$Size2 - d2$Overlap
    m = matrix(c(TP,FP,FN,TN),nrow=2)
    f = fisher.test(m, conf.level = 0.95)
    dataset = rbind(dataset,  data.table(ExternalList = extList, odds = f$estimate, pval = f$p.value, confmin = f$conf.int[1], confmax = f$conf.int[2]) )
  }
  
  dataset$odds = log2(dataset$odds)
  dataset$confmin = log2(dataset$confmin)
  dataset$confmax = log2(dataset$confmax)
  
  dataset[odds < -MAX_X_AXIS]$odds = -MAX_X_AXIS
  dataset[odds > MAX_X_AXIS]$odds = MAX_X_AXIS
  dataset[confmin < -MAX_X_AXIS]$confmin = -MAX_X_AXIS
  dataset[confmax > MAX_X_AXIS]$confmax = MAX_X_AXIS
  dataset$idx = seq(1,nrow(dataset)*3,3)
  return(dataset)
}


data1 = "lof_ch_genes"
data2 = "lof_possible_no_ch_genes"
dataset1 = direct_fisher(data,data1,data2)

data1 = "lof_high_ch_genes"
data2 = "lof_high_possible_no_ch_genes"
dataset2 = direct_fisher(data,data1,data2)

data1 = "rand_ch_genes"
data2 = "rand_possible_no_ch_genes"
dataset3 = direct_fisher(data,data1,data2)

f1 = dataset1[c(1,2,3,4,5),]
f2 = dataset2[c(1,2,3,4,5),]
f3 = dataset3[c(1,2,3,4,5),]

factor = 0.8

g2 = ggplot() +
  geom_rect(data = f1[idx %% 2 == 0], aes(xmin = -Inf, xmax = Inf, ymin = idx-1.5, ymax = idx+1.5), fill = "#d9d9d9", size = 1) +
  geom_vline(xintercept = 0, linetype = 2, size = 0.7) +
  geom_segment(data = f1, aes(x = confmin, xend = confmax, y = idx+factor, yend = idx+factor ), color = "black", size = 3.5 ) + 
  geom_segment(data = f1, aes(x = confmin, xend = confmax, y = idx+factor, yend = idx+factor ), color = colors[1], size = 2.5 ) + 
  geom_point(data = f1, aes(x = odds, y = idx+factor ), fill = colors[1], size = 4, shape = 21) + 
  geom_segment(data = f2, aes(x = confmin, xend = confmax, y = idx, yend = idx ), color = "black", size = 3.5 ) + 
  geom_segment(data = f2, aes(x = confmin, xend = confmax, y = idx, yend = idx ), color = colors[2], size = 2.5 ) + 
  geom_point(data = f2, aes(x = odds, y = idx), fill = colors[2], size = 4, shape = 21) +
  geom_segment(data = f3, aes(x = confmin, xend = confmax, y = idx-factor, yend = idx-factor ), color = "black", size = 3.5 ) + 
  geom_segment(data = f3, aes(x = confmin, xend = confmax, y = idx-factor, yend = idx-factor ), color = colors[3], size = 2.5 ) +
  geom_point(data = f3, aes(x = odds, y = idx-factor), fill = colors[3], size = 4, shape = 21) + 
  xlab("log2(Odds ratio)") +
  ylab("") +
  labs(tag = "b") +
  scale_x_continuous(limits = c(-MAX_X_AXIS,MAX_X_AXIS)) +
  scale_y_continuous(breaks = seq(1,nrow(f1)*3,3), labels = f1$ExternalList) +
  theme_minimal() +
  theme(text = element_text(size=12), plot.tag = element_text(margin = margin(b = -10, r = -30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

#############
# Panel c
#############

d1 = dataset1[c(6,7,8),]
d2 = dataset2[c(6,7,8),]
d3 = dataset3[c(6,7,8),]

g3 = ggplot() +
  geom_rect(data = d1[idx == 19], aes(xmin = -Inf, xmax = Inf, ymin = idx-1.5, ymax = idx+1.5), fill = "#d9d9d9", size = 1) +
  geom_vline(xintercept = 0, linetype = 2, size = 0.7) +
  geom_segment(data = d1, aes(x = confmin, xend = confmax, y = idx+factor, yend = idx+factor ), color = "black", size = 3.5 ) + 
  geom_segment(data = d1, aes(x = confmin, xend = confmax, y = idx+factor, yend = idx+factor ), color = colors[1], size = 2.5 ) +
  geom_point(data = d1, aes(x = odds, y = idx+factor ), fill = colors[1], size = 4, shape = 21) + 
  geom_segment(data = d2, aes(x = confmin, xend = confmax, y = idx, yend = idx ), color = "black", size = 3.5 ) +
  geom_segment(data = d2, aes(x = confmin, xend = confmax, y = idx, yend = idx ), color = colors[2], size = 2.5 ) + 
  geom_point(data = d2, aes(x = odds, y = idx), fill = colors[2], size = 4, shape = 21) +
  geom_segment(data = d3, aes(x = confmin, xend = confmax, y = idx-factor, yend = idx-factor ), color = "black", size = 3.5 ) + 
  geom_segment(data = d3, aes(x = confmin, xend = confmax, y = idx-factor, yend = idx-factor ), color = colors[3], size = 2.5 ) + 
  geom_point(data = d3, aes(x = odds, y = idx-factor), fill = colors[3], size = 4, shape = 21)  +
  labs(tag = "c") +
  xlab("log2(Odds ratio)") +
  ylab("") +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_x_continuous(limits = c(-2,2)) +
  scale_y_continuous(breaks = c(16,19,22), limits = c(15,23), labels = d1$ExternalList) +
  theme_minimal() +
  theme(text = element_text(size=12), plot.tag = element_text(margin = margin(b = -10, r = -30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))


#############
# Panel d
#############

dataLoF$LoF = dataLoF$numCH/dataLoF$expectedCH
dataSyn$Syn = dataSyn$numCH/dataSyn$expectedCH
dataMis$Mis = dataMis$numCH/dataMis$expectedCH

m1 = merge(dataLoF[,.(gene,LoF)],dataSyn[,.(gene,Syn)], by = "gene")
mergedData = merge(m1, dataMis[,.(gene,Mis)], by = "gene")
meltedData = melt(mergedData)
meltedData = meltedData[!is.na(value)]

t = wilcox.test(meltedData[variable == "LoF"][!is.na(value)]$value,meltedData[variable == "Mis"][!is.na(value)]$value)
t$p.value
t = wilcox.test(meltedData[variable == "Syn"][!is.na(value)]$value,meltedData[variable == "Mis"][!is.na(value)]$value)
t$p.value

medianDT = melt(data.table(LoF = round(median(meltedData[variable == "LoF"]$value),1), Mis = median(meltedData[variable == "Mis"]$value), Syn = round(median(meltedData[variable == "Syn"]$value),1)) )
medianDT$value2 = c("0.0","0.8","1.4") # forcing 0 to be 0.0

g4 = ggplot(meltedData,aes(x = variable, y = value, fill = variable)) + 
  geom_boxplot(width = 0.4, size = 0.5, outlier.size = 0.2) +
  geom_text(data = medianDT, aes(x = variable, y = value, label = value2 ), size = 3, hjust = -1.5, fontface = "bold") +
  geom_signif(comparisons = list(c("LoF", "Mis"),c("Mis", "Syn")), textsize = 3, step_increase = 0.04, tip_length = 0.02, annotation = c("p<5e-324","p=1e-170")) + 
  labs(tag = "d") +
  xlab("Variant type") +
  ylab("Observed / Expected compound hets.") +
  scale_x_discrete( limits = c("LoF","Mis","Syn"), labels = c("LoF","Missense","Synonymous")) +
  scale_fill_manual(values = c("#737373","#f0f0f0","#bdbdbd")) +
  theme_minimal() +
  theme(text = element_text(size = 12), plot.title = element_text(hjust = 0.5), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position = "none",
        panel.background = element_rect(colour = "black", fill=NA, size = 1)  ) # aspect.ratio = 1


# # ## DRAW LEGEND
legendDT = data.table(cols = colors, data = labels, x = c(1,2,3), y = c(1,2,3))
g5 = ggplot(legendDT, aes(x = x, y = y, fill = data)) +
  geom_point(size = 3.5, shape = 21, color = "black") +
  scale_fill_manual(values = legendDT$cols, labels = labels) +
  guides(fill = guide_legend(keywidth = 0.9,  keyheight = 0.9)) + 
  theme_minimal() +
  theme(text = element_text(size = 14), legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"), legend.spacing.y = unit(0, "mm"), legend.spacing.x = unit(0, "mm"), legend.key=element_blank(), legend.text = element_text(margin = margin(t = 0)))

legend <- cowplot::get_legend(g5)


design <- "
1113333
1113333
1113333
1113333
1113333
1113333
1113333
1113333
2223333
2224444
2224444
2224444
2224444
2224444"

g1 + inset_element(legend, left = 0.6, bottom = 1, right = 0.0, top = 0.5, on_top = FALSE, align_to = "panel") + g4 + g2 + inset_element(legend, left = 0.5, bottom = 1, right = 1, top = 0, on_top = FALSE, align_to = "panel") + g3 + plot_layout(design = design)


dev.off()
