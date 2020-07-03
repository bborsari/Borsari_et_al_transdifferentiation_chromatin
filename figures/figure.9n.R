.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")



#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(ggpubr)
library(scales)
library(xtable)



#********
# BEGIN *
#********

# 1. set working directory
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


# 2. read expression matrix
expression.m <- read.table("expression/QN.merged/expression.matrix.tsv",
                           h=T, sep="\t", stringsAsFactors = F)

# 3. keep expression values for 0h and 168h
expression.m <- expression.m[, c("H000", "H168")]


# 4. compute fold-change between 0h and 168h
expression.m$FC <- expression.m$H168 - expression.m$H000
expression.m$gene_id <- rownames(expression.m)


# 5. retrieve sets of up-regulated genes
metadata <- read.table("H3K4me3/QN.merged/metadata.tsv", h=T, sep="\t",
                       stringsAsFactors = F)
metadata <- metadata[metadata$final_class == "upregulation", ]
metadata$gene_id <- rownames(metadata)


# 6. subset expression matrix for up-regulated genes
expression.m <- merge(expression.m, metadata[, c("gene_id", "final_class")],
                      by = "gene_id")


# 7. read clusters df
clusters <- read.table("../Hi-C/cluster.1.2.3.txt", h=F, sep="\t",
                       stringsAsFactors = F)
colnames(clusters) <- c("gene_id", "cluster")


# 8. add to expression matrix clusters info
expression.m <- merge(expression.m, clusters, by = "gene_id")
colnames(expression.m)[2:4] <- c("0 h", "168 h", "fold-change")


# 9. prepare df for plot
df.plot <- melt(expression.m, id.vars = c("gene_id", "final_class", "cluster"))


# 10. define color palette
palette <- c("#F56E47", "#97BF04", "#772B59")
names(palette) <- c("cluster1", "cluster2", "cluster3")


# 11. make plots
x <- c("0 h", "168 h", "fold-change")

lop <- list()


for (i in 1:3) {
  
  lop[[i]] <- 
    ggplot(df.plot[df.plot$variable == x[i], ],
           aes(x=cluster, y=value, fill=cluster)) +
    geom_violin(alpha=.4, colour="white", scale = "width") +
    geom_boxplot(width=0.25, alpha = .6, outlier.shape = NA) +
    stat_compare_means(comparisons = list(c("cluster1", "cluster2"),
                                          c("cluster2", "cluster3"),
                                          c("cluster1", "cluster3"))) +
    scale_x_discrete(labels = c("1", "2", "3")) +
    ylim(0, 20) +
    scale_fill_manual(values = palette) +
    theme_bw() +
    theme(axis.title.x = element_text(size=12),
          axis.title.y = element_text(size=12),
          plot.title = element_text(size = 15, hjust = .5),
          axis.text.y = element_text(size =12),
          axis.text.x = element_text(size=12),
          strip.text.x = element_text(size=12),
          strip.background.x = element_blank(),
          panel.border = element_rect(color="black"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.title = element_blank()) +
    guides(fill=F) +
    labs(title = x[i])
  

  if ( i == 1 ) { lop[[i]] <- lop[[i]] + ylab("log2(TPM+1)") } 
  else {lop[[i]] <- lop[[i]] + ylab("")}

  if ( i == 2 ) { lop[[i]] <- lop[[i]] + xlab("clusters") }
  else {lop[[i]] <- lop[[i]] + xlab("")}
  
  
  
}



# 12. save plots
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.9n.pdf",
    height = 4, width = 8)
plot_grid(plotlist = lop, nrow=1, align = "h")
dev.off()

