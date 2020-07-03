.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(pheatmap)
library(reshape2)


#************
# FUNCTIONS *
#************

palette <- c("no_peak" = "black",
             "mark_ant" = "#a1d99b",
             "no_diff" = "#41ab5d",
             "exp_ant" = "#006d2c",
             "stable" = "#fec44f",
             "not_correlated" = "#ec7014",
             "negatively_correlated" = "#993404")


hours <- c("1" = 0, "2" = 3, "3" = 6, "4" = 9,
           "5" = 12, "6" = 18, "7" = 24,
           "8" = 36, "9" = 48, "10" = 72,
           "11" = 120, "12" = 168)

#********
# BEGIN *
#********

expression.m <- read.table("H3K4me3/QN.merged/expression.matrix.tsv", h=T, sep="\t")
expression.m$min <- apply(expression.m, 1, min)
expression.m$max <- apply(expression.m, 1, max)
expression.m <- expression.m[, c("min", "max")]


cluster1 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/GO.analysis/upregulated.1.genes.txt",
                       h=F, sep="\t", stringsAsFactors = F)
cluster3 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/GO.analysis/upregulated.3.genes.txt",
                       h=F, sep="\t", stringsAsFactors = F)
cluster2.lowest <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/GO.analysis/upregulated.2.lowest.txt",
                              h=F, sep="\t", stringsAsFactors = F)
cluster2.highest <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/GO.analysis/upregulated.2.highest.txt",
                              h=F, sep="\t", stringsAsFactors = F)

cluster1 <- cluster1$V1
cluster3 <- cluster3$V1
cluster2.lowest <- cluster2.lowest$V1
cluster2.highest <- cluster2.highest$V1

expression.m <- expression.m[rownames(expression.m) %in% 
                               c(cluster1, cluster3, cluster2.highest, cluster2.lowest),]

groups <- c()

for ( i in 1:nrow(expression.m) ) {
  
  if (rownames(expression.m)[i] %in% cluster1) {
    
    groups <- c(groups, "cluster1")
    
  } else if (rownames(expression.m)[i] %in% cluster3) {
    
    groups <- c(groups, "cluster3")
    
  } else if (rownames(expression.m)[i] %in% cluster2.lowest) {
    
    groups <- c(groups, "cluster2.lowest")
    
  } else {
    
    groups <- c(groups, "cluster2.highest")
    
  }
  
}


expression.m$group <- groups
expression.m$FC <- expression.m$max - expression.m$min
expression.m.melt <- melt(expression.m)
expression.m.melt$group <- factor(expression.m.melt$group,
                                  levels = c("cluster3",
                                             "cluster2.lowest",
                                             "cluster2.highest",
                                             "cluster1"))

expression.m$group <- factor(expression.m$group,
                             levels = c("cluster3",
                                        "cluster2.lowest",
                                        "cluster2.highest",
                                        "cluster1"))

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8d.1.pdf",
    height= 4, width =4)
ggplot(expression.m, aes(x=group, y=min, fill=group)) +
  stat_compare_means(comparisons = list(c("cluster3", "cluster2.lowest"),
                                        c("cluster2.lowest", "cluster2.highest"),
                                        c("cluster1", "cluster2.highest"),
                                        c("cluster3", "cluster2.highest"),
                                        c("cluster2.lowest", "cluster1"),
                                        c("cluster3", "cluster1")),
                     size = 3.5) + 
  geom_boxplot(width=0.5) +
  guides(fill=F) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        plot.title = element_text(size = 15, hjust = .5),
        axis.text.y = element_text(size =15),
        axis.text.x = element_text(size=12, angle = 30, vjust = .5),
        strip.text.x = element_text(size=15),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank()) +
  labs(title = "minimum expression") +
  ylab("log2(TPM+1)") +
  ylim(0, 25)
dev.off()


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8d.2.pdf",
    height= 4, width =4)
ggplot(expression.m, aes(x=group, y=max, fill=group)) +
  stat_compare_means(comparisons = list(c("cluster3", "cluster2.lowest"),
                                        c("cluster2.lowest", "cluster2.highest"),
                                        c("cluster1", "cluster2.highest"),
                                        c("cluster3", "cluster2.highest"),
                                        c("cluster2.lowest", "cluster1"),
                                        c("cluster3", "cluster1")),
                     size = 3.5) + 
  geom_boxplot(width=0.5) +
  guides(fill=F) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        plot.title = element_text(size = 15, hjust = .5),
        axis.text.y = element_text(size =15),
        axis.text.x = element_text(size=12, angle = 30, vjust = .5),
        strip.text.x = element_text(size=15),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank()) +
  labs(title = "maximum expression") +
  ylab("log2(TPM+1)") +
  ylim(0,25)
dev.off()





pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8d.3.pdf",
    height= 4, width =4)
ggplot(expression.m, aes(x=group, y=FC, fill=group)) +
  stat_compare_means(comparisons = list(c("cluster3", "cluster2.lowest"),
                                        c("cluster2.lowest", "cluster2.highest"),
                                        c("cluster1", "cluster2.highest"),
                                        c("cluster3", "cluster2.highest"),
                                        c("cluster2.lowest", "cluster1"),
                                        c("cluster3", "cluster1")),
                     size = 3.5) + 
  geom_boxplot(width=0.5) +
  guides(fill=F) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        plot.title = element_text(size = 15, hjust = .5),
        axis.text.y = element_text(size =15),
        axis.text.x = element_text(size=12, angle = 30, vjust = .5),
        strip.text.x = element_text(size=15),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank()) +
  labs(title = "fold-change") +
  ylab("log2(TPM+1)") 
dev.off()



tp <- read.table("H3K4me3/QN.merged/ant.del.analysis/increases/expression.25.50.75.100.tsv", 
                 h=T, sep="\t")
tp <- as.data.frame(t(apply(tp, 1, function(x){x <- hours[x]})))
colnames(tp) <- paste0(rep("perc_", 4), c(25,50,75,100))

tp.cluster1 <- tp[rownames(tp) %in% cluster1, ]
tp.cluster1$cluster <- "cluster1"

tp.cluster3 <- tp[rownames(tp) %in% cluster3, ]
tp.cluster3$cluster <- "cluster3"

tp.cluster2.highest <- tp[rownames(tp) %in% cluster2.highest, ]
tp.cluster2.highest$cluster <- "cluster2.highest"

tp.cluster2.lowest <- tp[rownames(tp) %in% cluster2.lowest, ]
tp.cluster2.lowest$cluster <- "cluster2.lowest"

tp <- as.data.frame(rbind(tp.cluster1,
                          tp.cluster2.lowest,
                          tp.cluster2.highest,
                          tp.cluster3))

tp$cluster <- factor(tp$cluster,
                     levels = c("cluster3",
                                "cluster2.lowest",
                                "cluster2.highest",
                                "cluster1"))

tp.melt <- melt(tp)


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8d.4.pdf",
    height = 4, width = 12)
ggplot(tp.melt, 
       aes(x=cluster, y=value)) +
  geom_boxplot(width=0.5) +
  guides(fill=F) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        plot.title = element_text(size = 15, hjust = .5),
        axis.text.y = element_text(size =15),
        axis.text.x = element_text(size=12, angle = 30, vjust = .5),
        strip.text.x = element_text(size=15),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank()) +
  ylab("time (hours)") +
  ylim(0, 168) +
  facet_grid(~variable)
dev.off()
