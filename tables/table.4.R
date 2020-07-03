.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")

palette <- c("positively_correlated" = "#DEA450",
             "no_peak" = "#9B461F",
             "negatively_correlated" = "#a095a0",
             "stable" = "#5d6f62",
             "not_correlated" = "#42354C")

palette2 <- c("down-regulated" = "#2d7f89", 
              "bending" = "#7acbd5",
              "up-regulated" = "#89372d",
              "peaking" = "#d5847a",
              "flat" = "#8f8f8f")



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



#************
# FUNCTIONS *
#************

plot.min.max.expression <- function(class) {
  
  cl1.expression.class <- cl1.expression[rownames(cl1.expression) %in% rownames(metadata[metadata$final_class==class, ]), ]
  cl2.expression.class <- cl2.expression[rownames(cl2.expression) %in% rownames(metadata[metadata$final_class==class, ]), ]
  cl3.expression.class <- cl3.expression[rownames(cl3.expression) %in% rownames(metadata[metadata$final_class==class, ]), ]
  
  cl1.expression.class$cluster <- "1"
  cl2.expression.class$cluster <- "2"
  cl3.expression.class$cluster <- "3"
  
  cl1.expression.class$minimum <- apply(cl1.expression.class[, 1:12], 1, min)
  cl2.expression.class$minimum <- apply(cl2.expression.class[, 1:12], 1, min)
  cl3.expression.class$minimum <- apply(cl3.expression.class[, 1:12], 1, min)
  
  cl1.expression.class$maximum <- apply(cl1.expression.class[, 1:12], 1, max)
  cl2.expression.class$maximum <- apply(cl2.expression.class[, 1:12], 1, max)
  cl3.expression.class$maximum <- apply(cl3.expression.class[, 1:12], 1, max)
  
  cl1.expression.class$fold_change <- cl1.expression.class$max - cl1.expression.class$min
  cl2.expression.class$fold_change <- cl2.expression.class$max - cl2.expression.class$min
  cl3.expression.class$fold_change <- cl3.expression.class$max - cl3.expression.class$min
  
  
  all.clusters.expression.class <- rbind(cl1.expression.class,
                                         cl2.expression.class,
                                         cl3.expression.class)
  
  all.clusters.expression.class.melt <- melt(all.clusters.expression.class)
  
  p <- ggplot(all.clusters.expression.class.melt[all.clusters.expression.class.melt$variable %in% c("minimum", "maximum", "fold_change"), ], 
              aes(x=cluster, y=value, fill=cluster)) + 
    geom_violin(alpha=.4, colour="white", width = 1) +
    geom_boxplot(width=0.25, alpha = .6, outlier.shape = NA) +
    guides(fill=F) +
    theme_bw() + 
    theme(axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          plot.title = element_text(size = 15, hjust = .5),
          axis.text.y = element_text(size =15),
          axis.text.x = element_text(size=15),
          strip.text.x = element_text(size=15),
          strip.background.x = element_blank(),
          panel.border = element_rect(color="black"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.title = element_blank()) +
    labs(title = class) +
    ylab("log2(TPM+1)") +
    xlab("clusters") +
    stat_compare_means(comparisons = list(c("1", "2"),
                                          c("1", "3"),
                                          c("2", "3")),
                       method = "wilcox.test",
                       size = 4.5) +
    facet_grid(~variable) +
    ylim(0, 21) +
    scale_fill_manual(values = c("#F56E47", "#97BF04", "#772B59"))
  
  return(p)
  
}




#********
# BEGIN *
#********

# 1. read dataframes

cl1 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.1.txt", h=F, sep="\t")
cl2 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.2.txt", h=F, sep="\t")
cl3 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.3.txt", h=F, sep="\t")
expression.matrix <- read.table("H3K4me3/QN.merged/expression.matrix.tsv", h=T, sep="\t")


colnames(cl1) <- "gene_id"
colnames(cl2) <- "gene_id"
colnames(cl3) <- "gene_id"


cl1.expression <- expression.matrix[rownames(expression.matrix) %in% cl1$gene_id, ]
cl2.expression <- expression.matrix[rownames(expression.matrix) %in% cl2$gene_id, ]
cl3.expression <- expression.matrix[rownames(expression.matrix) %in% cl3$gene_id, ]

classes <- c("up-regulated", "peaking","down-regulated", "bending")


# 2. read metadata file

metadata <- read.table("H3K4me3/QN.merged/metadata.tsv", h=T, sep="\t")
metadata$gene_id <- rownames(metadata)
metadata$final_class <- gsub("regulation", "-regulated", metadata$final_class)

lop <- list()

# 3. retrieve number of genes belonging to the different classes within each cluster

## cluster #1
cl1 <- merge(cl1, metadata[, c("gene_id", "final_class")], by = "gene_id")
cl1.v <- c("1", table(cl1$final_class))

## cluster #2
cl2 <- merge(cl2, metadata[, c("gene_id", "final_class")], by = "gene_id")
cl2.v <- c("2", table(cl2$final_class))

## cluster #3
cl3 <- merge(cl3, metadata[, c("gene_id", "final_class")], by = "gene_id")
cl3.v <- c("3", table(cl3$final_class))



# 4. generate plots (of % of classes per cluster)
all.clusters <- as.data.frame(rbind(cl1.v, cl2.v, cl3.v))
colnames(all.clusters)[1] <- "cluster"
all.clusters <- apply(all.clusters, 2, as.character)

for (x in 1:3) {
  
  z = sum(as.numeric(t(all.clusters[x, 2:5])))

  for (y in 2:5) {
    
    n = as.numeric(as.character(all.clusters[x,y]))
  
    all.clusters[x,y] <- paste0(n, " (", round((n/z)*100), "%)")
          
  }
  
}


print(xtable(all.clusters), include.rownames=FALSE)
