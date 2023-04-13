setwd('~/Dropbox/SHAPEIT5/Final_revision/Figures/Scripts/extended_figure7/')
#29-Sep-2022 Diogo Ribeiro @ UNIL
# CH results summary

library(data.table)
library(ggplot2)

jpeg("../../Extended_figures/Extended_figure7.jpeg", 1800, 1400, quality = 100, res=300)

data1 = fread("../../Source_data/CH/lof.nosingle.results.txt")
data2 = fread("../../Source_data/CH/lof.nosingle.cutoff0.99.results.txt")
data3 = fread("../../Source_data/CH/lof.nosingle.rand.results.txt")

labels = c("Full data","High confidence","Random phasing")
colors = c("#c6dbef","#2171b5","#fc9272")

#############
# Plot # CH events per category
#############

c1CH = sum(data1$numCH)
c2CH = sum(data2$numCH)
c3CH = sum(data3$numCH)
dt = data.table(labels = labels, CH = c(c1CH,c2CH,c3CH), cols = colors)

ggplot(dt, aes(x = labels, y = CH, fill = cols)) + 
  geom_bar(stat = "identity", size = 1, color = "black", width = 0.6) +
  geom_text(data = dt[c(1,2,3,4),], aes(label = CH), size = 4, vjust = -0.3, fontface = "bold") +
  geom_text(data = dt[c(5),], aes(label = CH), size = 4, vjust = 1.3, fontface = "bold") +
  theme_minimal() +
  scale_fill_identity(labels = dt$labels, guide = guide_legend(reverse = TRUE)) + 
  scale_x_discrete(limits = dt$labels) +
  ylab("# compound heterozygous genes") +
  xlab("Dataset") +
  theme(text = element_text(size = 14), plot.title = element_text(hjust = 0.5), 
        legend.position = "none",
        legend.box.background = element_rect(colour = "black"), legend.spacing.y = unit(0, "mm"), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1),  ) #aspect.ratio = 1 # axis.text.x = element_text(angle = 20, vjust=0.6)

dev.off()
