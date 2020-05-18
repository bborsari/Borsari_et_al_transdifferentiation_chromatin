.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(reshape2)
library(ggalluvial)
library(viridis)
library(ggridges)
library(ggpubr)




#********
# BEGIN *
#********

# 1. set working directory 
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


# 2. define correspondence between hours and time-points
hours <- c("1" = 0, "2" = 3, "3" = 6, "4" = 9,
           "5" = 12, "6" = 18, "7" = 24,
           "8" = 36, "9" = 48, "10" = 72,
           "11" = 120, "12" = 168)

# 3. read expression dynamics df
m <- read.table("H3K4me3/QN.merged/ant.del.analysis/decreases/expression.75.50.25.0.tsv",
                h=T, sep="\t")
colnames(m) <- c("75%", "50%", "25%", "0%")
m$gene_id <- rownames(m)


# 4. substitute time-points with hours
m[, 1:4] <- apply(m[, 1:4], 2, function(x){x <- hours[x]})


# 5. read clusters df
clusters <- read.table("../Hi-C/cluster.1.2.3.txt", h=F, sep="\t")
colnames(clusters) <- c("gene_id", "cluster")


# 6. subset clusters df for down-regulated genes
metadata <- read.table("H3K4me3/QN.merged/metadata.tsv", h=T, sep="\t",
                       stringsAsFactors = F)
clusters <- clusters[clusters$gene_id %in% 
                       rownames(metadata[metadata$final_class == "downregulation", ]), ]


# 7. merge m and clusters
clusters <- merge(clusters, m, by="gene_id")
rownames(clusters) <- clusters$gene_id


# 8. prepare df for plot
clusters.melt <- melt(clusters[, 2:6], id.vars = "cluster")
clusters.melt$value <- factor(clusters.melt$value, levels = hours)


# 9. define color palette
palette <- c("cluster1" = "#F56E47",
             "cluster2" = "#97BF04", 
             "cluster3" = "#772B59")

# 9. make plot
lop <- list()
x <- c("75%", "50%", "25%", "0%")

for (i in 1:4) {
  
  lop[[i]] <- ggplot(clusters.melt[clusters.melt$variable == x[i], ], 
                     aes(x=value, fill=cluster)) +
    geom_bar() +
    facet_grid(cluster~variable, scales = "free_y") +
    scale_x_discrete(drop=FALSE) +
    scale_fill_manual(values = palette) +
    guides(fill=F) +
    theme_bw() +
    theme(axis.title.x = element_text(size=8, hjust = .5),
          axis.title.y = element_text(size=8, vjust = .5),
          axis.text.x = element_text(size=10, angle = 30, vjust=.5),
          axis.text.y = element_text(size=13),
          strip.text.x = element_text(size=15),
          axis.ticks.x = element_line(color="black"),
          panel.border = element_rect(color="black"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_blank(),
          legend.title = element_blank(),
          legend.text=element_text(size=15),
          legend.position = "bottom",
          strip.background = element_blank(),
          strip.text = element_blank(),
          plot.title = element_blank()) +
    xlab("time (hours)") 
  
  if (i == 1) {lop[[i]] <- lop[[i]] + ylab("number of genes")} 
  else {lop[[i]] <- lop[[i]] + ylab("")}
  
}




pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.9j.pdf", 
    height=8, width=12)
plot_grid(plotlist = lop, nrow = 1, ncol=4)
dev.off()
