.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")



#************
# LIBRARIES *
#************


library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(circlize)
library(seriation)
library(ggplot2)
library(ggpubr)
library(cowplot)


#************
# FUNCTIONS *
#************

function1 <- function(mark) {
  
  # 1. read mark matrix of 12448 genes after QN merged normalization
  mark.matrix <- read.table(paste0(mark, "/QN.merged/", mark, ".matrix.after.QN.merged.tsv"),
                            h=T, sep="\t")
  
  
  # 2. retrieve set of marked genes
  marked.genes <- read.table(paste0(mark, "/QN.merged/all.genes.intersecting.peaks.tsv"),
                             h=F, sep="\t", stringsAsFactors = F)
  marked.genes <- marked.genes$V1
  
  
  # 3. retrieve set of genes w/ variable profile for the mark
  mark.sig.genes <- read.table(paste0(mark, "/QN.merged/", mark, ".QN.merged.maSigPro.out.tsv"), 
                               h=T, sep="\t")
  mark.sig.genes$gene_id <- rownames(mark.sig.genes)
  mark.sig.genes <- mark.sig.genes %>% separate(gene_id, c("gene", "id"), "\\.")
  mark.sig.genes <- mark.sig.genes$gene

  
  # 4. retrieve differentially marked genes (aka w/ peak of the mark and variable profile)
  mark.sig.genes <- intersect(mark.sig.genes, marked.genes)
  
  print(mark)
  print(length(mark.sig.genes))
  
  
  # 5. compute fold-change
  x <- mark.matrix[rownames(mark.matrix) %in% mark.sig.genes, ]
  x$FC <- apply(x, 1, function(x){max(x) - min(x)})
  
  
  # 6. add mark and gene_id
  x$mark <- mark
  x$gene_id <- rownames(x)
  
  
  # 7. add group of genes
  x <- merge(x[, c("gene_id", "FC", "mark")], m, all.x = T, by = "gene_id")
  
  
  # 6. return mark matrix scaled for differentially marked genes
  return(x)
  
  
}


#********
# BEGIN *
#********


# 1. set working directory
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/")


# 2. read expression matrix of 12448 genes
expression.matrix <- read.table("expression/QN.merged/selected.genes.rep.2.3.after.QN.merged.tsv",
                                h=T, sep="\t")


# 3. retrieve set of not expressed genes
not.expressed.genes <- read.table("expression/silent.genes.txt", h=F, sep="\t", 
                                  stringsAsFactors = F)
not.expressed.genes <- not.expressed.genes$V1


# 4. retrieve set of DE genes
DE.genes <- read.table("expression/QN.merged/expression.matrix.tsv", h=T, sep="\t")
DE.genes <- rownames(DE.genes)


# 5. retrieve set of stably expressed genes
stable.genes <- setdiff(rownames(expression.matrix), c(not.expressed.genes, DE.genes))


# 6. create a df with info of not expressed, stably expressed and DE genes
m <- data.frame(gene_id = c(not.expressed.genes, stable.genes, DE.genes),
                group = c(rep("not_expressed", length(not.expressed.genes)),
                          rep("stable", length(stable.genes)),
                          rep("DE_genes", length(DE.genes))))


# 7. the marks we're analyzing
marks <- c("H3K27ac", "H3K9ac", "H4K20me1", "H3K4me3", "H3K4me1", "H3K36me3",
           "H3K4me2", "H3K9me3", "H3K27me3")


# 8. prepare df for plot
df.plot <- data.frame(stringsAsFactors = F)
for ( i in 1:9 ) {
  
  df.plot <- rbind(df.plot, function1(mark = marks[i]))
  
}
df.plot$group <- factor(df.plot$group, levels = c("DE_genes", "stable", "not_expressed"))


# 9. make plots
lop <- list()

for ( i in 1:9 ) {
  
  lop[[i]] <- ggplot(df.plot[df.plot$mark == marks[i], ], aes(x=group, y=FC)) +
    geom_violin(aes(fill = group), color="black") +
    geom_boxplot(width = 0.3, outlier.shape = NA) +
    scale_fill_manual(values = c("DE_genes" = "#f768a1", 
                                 "stable" = "#bdbdbd", 
                                 "not_expressed" = "#737373")) +
    stat_compare_means( comparisons = list(c("DE_genes", "stable"),
                                           c("DE_genes", "not_expressed"),
                                           c("stable", "not_expressed")),
                        size = 3.5) +
    labs(title = marks[i]) +
    theme_bw() +  
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=10),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          panel.border = element_rect(color="black"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          plot.title = element_text(size=13, hjust = .5),
          axis.line = element_line(colour = "black")) +
    expand_limits(y = 0) +
    guides(fill=F) +
    scale_x_discrete(labels = c("DE", "stably e.", "not e.")) 
  
}



pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.9m.pdf",
    width = 15, height = 7)
plot_grid(plotlist = lop, nrow=2, align="hv")
dev.off()

