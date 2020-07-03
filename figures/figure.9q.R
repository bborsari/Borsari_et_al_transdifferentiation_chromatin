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
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/ENCODE.data")


# 2. read cell lines
cell.lines <- c("MCF-7", "HepG2", "A549", "GM12878", "K562",  
                "B.cells", "CD14.positive.monocytes")


# 3. read, for each cell line, the quantification df
x <- data.frame(stringsAsFactors = F)
for (i in 1:7) {
  
  # 3.1. read rep1
  rep1 <- read.table(paste0(cell.lines[i], 
                            "/quantifications/", 
                            cell.lines[i], 
                            ".1.tsv"), h=T, sep="\t")
  
  # 3.2. read rep2
  rep2 <- read.table(paste0(cell.lines[i], 
                            "/quantifications/", 
                            cell.lines[i], 
                            ".2.tsv"), h=T, sep="\t")
  
  # 3.3. merge replicates
  tmp <- merge(rep1[, c("gene_id", "TPM")], rep2[, c("gene_id", "TPM")],
               by = "gene_id")
  
  # 3.4. compute mean expression between replicates
  tmp$TPM <- apply(tmp[, 2:3], 1, mean)
  tmp$gene_id <- gsub("\\..*", "", tmp$gene_id)
  
  
  # 3.5. add cell line info
  tmp$cell_line <- cell.lines[i]
  x <- rbind(x, tmp[, c("gene_id", "TPM", "cell_line"), drop=F])
  
}


# 4. read clusters df
clusters <- read.table("../Hi-C/cluster.1.2.3.txt", h=F, sep="\t",
                       stringsAsFactors = F)
colnames(clusters) <- c("gene_id", "cluster")



# 5. add clusters info
x <- merge(x, clusters, by = "gene_id", all.x = T)
x <- x[complete.cases(x), ] # keep only the 8030 genes within the 3 clusters


# 6. retrieve sets of down-regulated genes
metadata <- read.table("../all.marks/H3K4me3/QN.merged/metadata.tsv", h=T, sep="\t",
                       stringsAsFactors = F)
metadata <- metadata[metadata$final_class == "downregulation", ]


# 7. subset x for down-regulated genes
x <- x[x$gene_id %in% rownames(metadata), ]


# 8. define color palette
palette <- c("#F56E47", "#97BF04", "#772B59")
names(palette) <- c("cluster1", "cluster2", "cluster3")


# 9. define titles for plots
titles <- c("MCF-7", "HepG2", "A549", "GM12878", "K562",  
            "B cells", "CD14+ monocytes")


# 10. make plots
lop <- list()
for ( i in 1:7 ) {
  
  lop[[i]] <- 
    ggplot(x[x$cell_line == cell.lines[i], ],
           aes(x=cluster, y=log2(TPM+1), fill=cluster)) +
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
    labs(title = titles[i]) 
  
  if ( i == 1 ) { lop[[i]] <- lop[[i]] + ylab("log2(TPM+1)") } 
  else {lop[[i]] <- lop[[i]] + ylab("")}
  
  if ( i == 4 ) { lop[[i]] <- lop[[i]] + xlab("clusters") } 
  else {lop[[i]] <- lop[[i]] + xlab("")}
  
  
}


# 11. save plots
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.9q.pdf",
    height = 4, width = 19)
plot_grid(plotlist = lop, nrow=1, ncol=7, align = "hv")
dev.off()

