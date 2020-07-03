.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")



palette <- c("H3K9ac" = "#af4d85",
             "H3K27ac" = "#630039",
             "H3K4me3" = "#d199b9",
             "H3K27me3" = "#1d2976",
             "H3K9me3" = "#a7add4",
             "H3K36me3" = "#7fbc41",
             "H4K20me1" = "#4c7027",
             "H3K4me1" = "#e5ab00",
             "H3K4me2" = "#a67c00",
             "expression" = "white")





#************
# LIBRARIES *
#************

library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(cowplot)
library(reshape2)
library(dplyr)
library(tidyr)


#************
# FUNCTIONS *
#************

my.function <- function(mark) {
  
  lop <- list()
  
  # 1. read mark matrix
  mark.matrix <- read.table(paste0(mark, "/QN.merged/", mark, ".matrix.after.QN.merged.tsv"),
                            h=T, sep="\t")
  
  # 2. check that rownames are sorted the same in expression and mark matrices
  stopifnot(identical(rownames(expression.matrix),
                      rownames(mark.matrix)))
  

  # 3. retrieve mark's specific significant genes 
  # and intersect with list of expression significant genes
  mark.sig.genes <- read.table(paste0(mark, "/QN.merged/", mark, ".QN.merged.maSigPro.out.tsv"), 
                               h=T, sep="\t")
  mark.sig.genes$gene_id <- rownames(mark.sig.genes)
  mark.sig.genes <- mark.sig.genes %>% separate(gene_id, c("gene", "id"), "\\.")
  rownames(mark.sig.genes) <- mark.sig.genes$gene
  mark.sig.genes$gene <- NULL
  mark.sig.genes$id <- NULL
  mark.sig.genes <- rownames(mark.sig.genes)
  expression.mark.sig.genes <- intersect(expression.sig.genes, mark.sig.genes)
  expression.mark.ns.genes <- setdiff(expression.sig.genes, expression.mark.sig.genes)
  
  # 4. filter for genes with peaks
  ## 4.1. genes w/ peaks
  mark.genes.peaks <- read.table(paste0(mark, "/QN.merged/all.genes.intersecting.peaks.tsv"),
                                 h=F, sep="\t", stringsAsFactors = F)
  mark.genes.peaks <- mark.genes.peaks$V1
  
  ## 4.2. filter sig and ns genes
  expression.mark.sig.genes <- intersect(expression.mark.sig.genes, mark.genes.peaks)
  expression.mark.ns.genes <- intersect(expression.mark.ns.genes, mark.genes.peaks)
  
  print(mark)
  print(length(expression.mark.sig.genes))
  print(length(expression.mark.ns.genes))
  
  # 5. get transposed matrix
  expression.matrix.f <- as.data.frame(t(expression.matrix))
  expression.matrix.scaled <- as.data.frame(scale(expression.matrix.f))
  
  mark.matrix <- as.data.frame(t(mark.matrix))
  mark.matrix.scaled <- as.data.frame(scale(mark.matrix))
  
  # 6. PCA with significantly variable
  # genes for both expression and mark
  stopifnot(identical(colnames(expression.matrix.scaled[, expression.mark.sig.genes]),
                      colnames(mark.matrix.scaled[, expression.mark.sig.genes])))
  
  pca.mark.sig <- prcomp(rbind(expression.matrix.scaled[, expression.mark.sig.genes],
                               mark.matrix.scaled[, expression.mark.sig.genes]),
                         center = F,
                         scale = F)
  
  df.pca.mark.sig <- as.data.frame(pca.mark.sig$x[, 1:4])
  
  summary.pca.mark.sig <- as.data.frame(summary(pca.mark.sig)$importance)
  
  df.pca.mark.sig$tp_type <- 
    c(paste0(rownames(df.pca.mark.sig)[1:12], "_expression"),
      paste0(rownames(df.pca.mark.sig)[1:12], "_mark"))
  
  df.pca.mark.sig$tp <- rep(rownames(df.pca.mark.sig)[1:12], 2)
  df.pca.mark.sig$type <- rep(c("expression", mark), each=12)
  df.pca.mark.sig$tp_numeric <- rep(c(0,3,6,9,12,18,24,36,48,72,120,168), 2)
  df.pca.mark.sig$group <- "differentially marked"
  
  
  lop[[1]] <- ggplot(df.pca.mark.sig, 
                     aes(x=PC1, y=PC2, fill = type)) +
    theme_bw() +
    facet_grid(~group) +
    theme(axis.title = element_text(size=12),
          axis.text = element_text(size=11),
          panel.border = element_rect(color="black"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          aspect.ratio = 1,
          strip.text.x = element_text(size=13),
          strip.background.x = element_blank(),
          legend.text = element_text(size=13),
          legend.title = element_text(size=13),
          plot.title = element_blank(),
          plot.margin = unit(rep(0, 4), "cm")) +    
    geom_point(aes(size=tp_numeric), shape=21, alpha=.85) +
    scale_fill_manual(values = palette) +
    coord_equal(ratio=1) +
    xlab(paste("PC1", round(summary.pca.mark.sig$PC1[2]*100, 1), "%")) +
    ylab(paste("PC2", round(summary.pca.mark.sig$PC2[2]*100, 1), "%")) +
    labs(title = mark) +
    guides(fill = F, size=F) 
  
  # 7. PCA with genes significantly
  # variable only for expression
  tmp <- as.data.frame(t(na.omit(t(mark.matrix.scaled[, expression.mark.ns.genes]))))
  expression.mark.ns.genes <- colnames(tmp)
  
  stopifnot(identical(colnames(expression.matrix.scaled[, expression.mark.ns.genes]),
                      colnames(mark.matrix.scaled[, expression.mark.ns.genes])))
  
  pca.mark.ns <- prcomp(rbind(expression.matrix.scaled[, expression.mark.ns.genes],
                              mark.matrix.scaled[, expression.mark.ns.genes]),
                           center = F,
                           scale = F)
  df.pca.mark.ns <- as.data.frame(pca.mark.ns$x[, 1:4])
  
  summary.pca.mark.ns <- as.data.frame(summary(pca.mark.ns)$importance)
  
  df.pca.mark.ns$tp_type <- 
    c(paste0(rownames(df.pca.mark.ns)[1:12], "_expression"),
      paste0(rownames(df.pca.mark.ns)[1:12], "_mark"))
  
  df.pca.mark.ns$tp <- rep(rownames(df.pca.mark.ns)[1:12], 2)
  df.pca.mark.ns$type <- rep(c("expression", mark), each=12)
  df.pca.mark.ns$tp_numeric <- rep(c(0,3,6,9,12,18,24,36,48,72,120,168), 2)
  df.pca.mark.ns$group <- "stable marking"
  
  lop[[2]] <- ggplot(df.pca.mark.ns, 
                        aes(x=PC1, y=PC2, fill = type)) +
    theme_bw() +
    facet_grid(~group) +
    theme(axis.title = element_text(size=12),
          axis.text = element_text(size=11),
          panel.border = element_rect(color="black"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          aspect.ratio = 1,
          strip.text.x = element_text(size=13),
          strip.background.x = element_blank(),
          legend.text = element_text(size=13),
          legend.title = element_text(size=13),
          plot.title = element_blank(),
          plot.margin = unit(rep(0, 4), "cm")) +
    geom_point(aes(size=tp_numeric), shape=21, alpha=.85) +
    scale_fill_manual(values = palette) +
    coord_equal(ratio=1) +
    xlab(paste("PC1", round(summary.pca.mark.ns$PC1[2]*100, 1), "%")) +
    ylab(paste("PC2", round(summary.pca.mark.ns$PC2[2]*100, 1), "%")) +
    guides(fill = F, size=F) 
  
  return(lop)
  
}



setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/")


#********
# BEGIN *
#********

# 1. read expression matrix
expression.matrix <- read.table("expression/QN.merged/selected.genes.rep.2.3.after.QN.merged.tsv",
                                h=T, sep="\t")

# 2. sort rownames of expression matrix
new.order <- sort(rownames(expression.matrix))
expression.matrix <- expression.matrix[new.order, ]

# 3. retrieve significantly variable genes
# for expression
# retrieve expression significant genes
expression.sig.genes <- read.table("expression/QN.merged/expression.QN.merged.maSigPro.out.tsv", 
                                   h=T, sep="\t")
expression.sig.genes <- rownames(expression.sig.genes)


# 5. make plots
main.lop <- list()
i=1

marks = c("H3K27ac", "H3K9ac", "H4K20me1",
          "H3K4me3", "H3K4me1", "H3K36me3",
          "H3K4me2", "H3K9me3", "H3K27me3")

for (my.mark in marks) {
  
  main.plot <- my.function(mark=my.mark) 
  main.lop[[i]]  <- plot_grid(main.plot[[1]], NULL, main.plot[[2]],
                              nrow=1, scale = c(.85, 0, .85), rel_widths = c(1, -0.25, 1),
                              labels = c(my.mark, "", ""))
  i <- i+1
}

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1g.pdf",
    height=11, width=21)
plot_grid(main.lop[[1]], NULL, main.lop[[2]], NULL, main.lop[[3]], 
          main.lop[[4]], NULL, main.lop[[5]], NULL, main.lop[[6]], 
          main.lop[[7]], NULL, main.lop[[8]], NULL, main.lop[[9]], 
          nrow=3, ncol=5, align = "hv", 
          rel_widths = c(1, -0.1, 1, -0.1, 1,
                         1, -0.1, 1, -0.1, 1,
                         1, -0.1, 1, -0.1, 1))
dev.off()
