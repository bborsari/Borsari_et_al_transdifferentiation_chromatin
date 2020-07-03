.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/nfs/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/Hi-C/")

#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(pheatmap)
library(reshape2)
library(dplyr)


#************
# FUNCTIONS *
#************

my.function1 <- function(f1.df, f1.title) {
  
  f1.df <- f1.df[, !(colnames(f1.df) %in% c("gene_id", "number_transitions"))]
  f1.df.melt <- melt(f1.df)
  
  
  f1.p <- ggplot(f1.df.melt, aes(x=cluster, y=value, fill=cluster)) +
    geom_boxplot(alpha = .7) +
    stat_compare_means(comparisons = list(c("cluster1", "cluster2"),
                                          c("cluster2", "cluster3"),
                                          c("cluster1", "cluster3")),
                       method = "wilcox.test",
                       paired = F) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=15),
          axis.text = element_text(size=15),
          plot.title = element_text(size=15, hjust = .5),
          panel.border = element_rect(color="black"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = c("#225ea8", "#7fcdbb", "#ff005a")) +
    labs(title = f1.title) +
    guides(fill=F) +
    scale_x_discrete(labels = c("1", "2", "3")) +
    ylab("compartment values")
  
  return(f1.p)
  
  
}



#********
# BEGIN *
#********

# 1. import complete datasets with A/B values for all time-points
# + number of transitions
all.tp <- read.table("all.time.points/complete.data.set.tsv", h=F, sep="\t")
matched.tp <- read.table("matching.time.points/complete.data.set.tsv", h=F, sep="\t")


# 2. update colnames
colnames(all.tp) <- c("gene_id", "0h", "9h",
                      "18h", "24h", "48h",
                      "72h", "96h", "120h",
                      "144h", "168h", "192h", "number_transitions")

colnames(matched.tp) <- c("gene_id", "0h", "9h",
                          "18h", "24h", "48h",
                          "72h", "120h",
                          "168h", "number_transitions")



# 3. import list of genes per cluster
clusters <- read.table("cluster.1.2.3.txt")

# 4. update colnames
colnames(clusters) <- c("gene_id", "cluster")

# 5. add to table of A/B values + transitions
# the info regarding the clusters
all.tp <- merge(all.tp, clusters, all.x = T, by = "gene_id")
matched.tp <- merge(matched.tp, clusters, all.x = T, by = "gene_id")


# 6. % of genes with at least 1 transition
lop <- list()
lop[[1]] <- my.function1(f1.df = all.tp, f1.title = "all time-points")
lop[[2]] <- my.function1(f1.df = matched.tp, f1.title = "matched time-points")

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8f.pdf",
    width = 7, height=4)
plot_grid(plotlist = lop, nrow=1, ncol=2, align="h")
dev.off()




