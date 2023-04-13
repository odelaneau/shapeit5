setwd('~/Dropbox/SHAPEIT5/Final_revision/Figures/Scripts/extended_figure9/')
#29-Sep-2022 Diogo Ribeiro @ UNIL
# CH results summary

library(data.table)
library(ggplot2)
library(patchwork)
library(cowplot)

jpeg("../../Extended_figures/Extended_figure9.jpeg", 3000, 2000, quality = 100, res=300)

dataLoF = fread("../../Source_data/CH/lof.nosingle.results.txt")
dataSyn = fread("../../Source_data/CH/synonymous.nosingle.results.txt")
dataMis = fread("../../Source_data/CH/missense.nosingle.results.txt")
data = fread("../../Source_data/CH/lof.enrichment.lof_mis_syn.txt", stringsAsFactors = FALSE, header = TRUE, sep="\t")

#############
# Panel a
#############

## Basic plot
nTotalGenes = nrow(dataSyn)
nLoFCH = nrow(dataLoF[numCH > 0])
nLoFPossible = nrow(dataLoF[numCH == 0][num2mut > 0])
nLoFOther = nTotalGenes - nLoFPossible - nLoFCH
nSynCH = nrow(dataSyn[numCH > 0])
nSynPossible = nrow(dataSyn[numCH == 0][num2mut > 0])
nSynOther = nTotalGenes - nSynPossible - nSynCH
nMisCH = nrow(dataMis[numCH > 0])
nMisPossible = nrow(dataMis[numCH == 0][num2mut > 0])
nMisOther = nTotalGenes - nMisPossible - nMisCH

basic = data.table(labels = c("LoF","Syn","Mis"), CH = c(nLoFCH,nSynCH,nMisCH))
g1 = ggplot(basic, aes(x = labels, y = CH, fill = labels)) +
  geom_bar(stat = "identity", size = 0.7, color = "black", width = 0.7) +
  geom_text(aes(label = CH), size = 3, vjust = -1, fontface = "bold") +
  labs(tag = "a") +
  ylab("# compound het. genes") +
  xlab("Variant type") +
  theme_minimal() +
  scale_x_discrete(limits = c("LoF","Mis","Syn"),labels = c("LoF","Missense","Synonymous")) +
  scale_fill_manual(values = c("#737373","#bdbdbd","#f0f0f0")) +
  ylim(c(0,max(basic$CH) + 1000)) +
  theme(text = element_text(size = 12), plot.title = element_text(hjust = 0.5),
          legend.position = "none", plot.tag = element_text(margin = margin(b = -10, r = -10)),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          panel.background = element_rect(colour = "black", fill=NA, size = 1) )


#############
# Panel b
#############

dt = data.table(group = c("LoF","LoF","LoF","Syn","Syn","Syn","Miss","Miss","Miss"), 
                Gene = c("3: ch","2: possible","1: other","3: ch","2: possible","1: other","3: ch","2: possible","1: other"), 
                values = c(nLoFCH,nLoFPossible,nLoFOther,nSynCH,nSynPossible,nSynOther,nMisCH,nMisPossible,nMisOther))

g2 = ggplot(dt,aes(x = group, fill = Gene, y = values)) + 
  geom_bar(stat = "identity", size = 0.7, color = "black", width = 0.6, alpha = 0.7) +
  labs(tag = "b") +
  scale_fill_manual(values = c("#66c2a5","#fc8d62","#8da0cb"), labels = c("Not assessed","Not Compound","Compound")) +
  scale_x_discrete(limits = c("Syn","Miss","LoF"), labels = c("Synonymous","Missense","LoF")) + 
  ylab("# Genes") +
  xlab("Variant type") +
  theme_minimal() +
  theme(text = element_text(size = 12), plot.title = element_text(hjust = 0.5), 
        plot.tag = element_text(margin = margin(b = -10, r = -10)),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill=NA, size = 1),
        legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=7), #change legend title font size
        legend.text = element_text(size=6), legend.position = "top", legend.margin=margin(b=-10))


#############
# Panel c
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

data1 = "mis_ch_genes"
data2 = "mis_possible_no_ch_genes"
dataset2 = direct_fisher(data,data1,data2)

data1 = "syn_ch_genes"
data2 = "syn_possible_no_ch_genes"
dataset3 = direct_fisher(data,data1,data2)

labels = c("LoF","Missense","Synonymous")
colors = c("#737373","#bdbdbd","#f0f0f0")

f1 = dataset1[c(1,2,3,4,5),]
f2 = dataset2[c(1,2,3,4,5),]
f3 = dataset3[c(1,2,3,4,5),]

factor = 0.8

g3 = ggplot() +
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
  labs(tag = "c") +
  scale_x_continuous(limits = c(-MAX_X_AXIS,MAX_X_AXIS)) +
  scale_y_continuous(breaks = seq(1,nrow(f1)*3,3), labels = f1$ExternalList) +
  theme_minimal() +
  theme(text = element_text(size=12), plot.tag = element_text(margin = margin(b = -10, r = -30)),
        # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

g3

#############
# Panel d
#############

d1 = dataset1[c(6,7,8),]
d2 = dataset2[c(6,7,8),]
d3 = dataset3[c(6,7,8),]

g4 = ggplot() +
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
  labs(tag = "d") +
  xlab("log2(Odds ratio)") +
  ylab("") +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_x_continuous(limits = c(-2,2)) +
  scale_y_continuous(breaks = c(16,19,22), limits = c(15,23), labels = d1$ExternalList) +
  theme_minimal() +
  theme(text = element_text(size=12), plot.tag = element_text(margin = margin(b = -10, r = -30)),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

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
2224444
"

# # ## DRAW LEGEND
# png("/home/dribeiro/legend_beagle.png",width = 600, height = 600)
legendDT = data.table(cols = colors, data = labels, x = c(1,2,3), y = c(1,2,3))
g5 = ggplot(legendDT, aes(x = x, y = y, fill = data)) +
  geom_point(size = 3.5, shape = 21, color = "black") +
  scale_fill_manual(values = legendDT$cols, labels = labels) +
  guides(fill = guide_legend(keywidth = 0.8,  keyheight = 0.8)) + 
  theme_minimal() +
  theme(text = element_text(size = 11), legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"), legend.spacing.y = unit(0, "mm"), legend.spacing.x = unit(0, "mm"), legend.key=element_blank(), legend.text = element_text(margin = margin(t = 0)))

legend <- cowplot::get_legend(g5)

g1 + inset_element(legend, left = 0.35, bottom = 1, right = 0.0, top = 0.5, on_top = FALSE, align_to = "panel") + g2 + g3 + inset_element(legend, left = 0.5, bottom = 1, right = 1, top = 0, on_top = FALSE, align_to = "panel") + g4 + plot_layout(design = design)


dev.off()

