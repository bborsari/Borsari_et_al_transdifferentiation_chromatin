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
all.clusters[, 2:ncol(all.clusters)] <- apply(all.clusters[, 2:ncol(all.clusters)], 2, as.numeric)
all.clusters.melt <- melt(all.clusters)
all.clusters.melt$variable <- factor(all.clusters.melt$variable, levels = c("up-regulated", 
                                                                            "peaking",
                                                                            "bending",
                                                                            "down-regulated"))

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3o1.pdf", width=4.5, height=3)
ggplot(all.clusters.melt, aes(x=cluster, y=value, fill=variable)) +
  geom_bar(stat="identity", width=0.7, position = "fill", colour="black") +
  theme_bw() + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size=15),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size=15),
        axis.text.x =element_text(angle=45, hjust=1, size=15),
        strip.text.x = element_text(size=15),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=palette2) + 
  scale_y_continuous(labels = percent_format()) +
  coord_flip() +
  ylab("% of genes") +
  xlab("clusters")
dev.off()


# 5. generate plots of min / max / FC per cluster


for (j in 1:length(classes)) {
  
  lop[[j]] <- plot.min.max.expression(class=classes[j])
  
}


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3o2.pdf", 
    width=6, height=4)
lop[[1]] <- lop[[1]] + theme(plot.title = element_blank())
lop[[1]]
dev.off()



pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3o3.pdf", 
    width=18, height=4.5)
plot_grid(plotlist = lop[2:4], nrow=1, ncol=3)
dev.off()




# # 4. multiple fisher tests between all groups and all classes
# pvals <- c()
# type <- c()
# 
# for (j in 2:(ncol(all.clusters))) {
#   
#   for (i in 1:3) {
#     
#     k = i+1
#     while (k <= 4) {
#       
#       aa = all.clusters[i, j]
#       bb = sum(all.clusters[i, 2:5]) - aa
#       cc = all.clusters[k, j]
#       dd = sum(all.clusters[k, 2:5]) - cc
#       
#       M <- matrix(c(aa, bb, cc, dd), byrow = T, nrow = 2)
#       
#       
#       type <- c(type, paste(colnames(all.clusters[j]), 
#                             as.character(all.clusters[i, 1]),
#                             as.character(all.clusters[k, 1])))
#       pvals <- c(pvals, fisher.test(M)$p.value)
#       
#       #print ( c(aa, bb, cc, dd))
#       k <- k+1
#       
#     }
#     
#     
#   }
#   
#   
#   
# }
# 
# 
# my.df <- data.frame(pvals = pvals, type = type)
# View(my.df[my.df$pvals < 0.05, ])
# 
# all.clusters[2, 2] / sum(all.clusters[2, 2:5]) * 100
# all.clusters[1, 2] / sum(all.clusters[1, 2:5]) * 100
# all.clusters[3, 2] / sum(all.clusters[3, 2:5]) * 100
# all.clusters[4, 2] / sum(all.clusters[4, 2:5]) * 100
# 
# all.clusters[1, 3] / sum(all.clusters[1, 2:5]) * 100
# all.clusters[2, 3] / sum(all.clusters[2, 2:5]) * 100
# all.clusters[3, 3] / sum(all.clusters[3, 2:5]) * 100
# all.clusters[4, 3] / sum(all.clusters[4, 2:5]) * 100
# 
# all.clusters[1, 4] / sum(all.clusters[1, 2:5]) * 100
# all.clusters[2, 4] / sum(all.clusters[2, 2:5]) * 100
# all.clusters[3, 4] / sum(all.clusters[3, 2:5]) * 100
# all.clusters[4, 4] / sum(all.clusters[4, 2:5]) * 100
# 
# all.clusters[1, 5] / sum(all.clusters[1, 2:5]) * 100
# all.clusters[2, 5] / sum(all.clusters[2, 2:5]) * 100
# all.clusters[3, 5] / sum(all.clusters[3, 2:5]) * 100
# all.clusters[4, 5] / sum(all.clusters[4, 2:5]) * 100
# 
# 
