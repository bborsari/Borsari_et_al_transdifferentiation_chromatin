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


my.function <- function(mark, peaks = F) {
  
  # 1. retrieve genes significantly variable for the mark
  mark.sig.genes <- read.table(paste0(mark, "/QN.merged/", mark, ".QN.merged.maSigPro.out.tsv"), 
                               stringsAsFactors = F, h=T,
                               sep="\t")
  rownames(mark.sig.genes) <- gsub("\\..*", "", rownames(mark.sig.genes))
  mark.sig.genes <- rownames(mark.sig.genes)
  
  if (peaks) {
    
    # 2. intersect these genes with genes that have a peak in at least one time point
    mark.genes.peaks <- read.table(paste0(mark, "/QN.merged/all.genes.intersecting.peaks.tsv"), 
                                   stringsAsFactors = F, h=F,
                                   sep="\t")
    
    mark.genes.peaks <- mark.genes.peaks$V1
    mark.sig.genes <- intersect(mark.sig.genes, mark.genes.peaks)
    
  }
  
  
  # 3. retrieve genes variable for expression
  variable.genes.expression <- read.table(paste0(mark, "/QN.merged/", mark, ".6.groups.tsv"),
                                          h=T, sep="\t")
  variable.genes.expression <- rownames(variable.genes.expression)
  
  
  # 4. retrieve silent genes
  silent.genes <- read.table("expression/silent.genes.txt", h=F, sep="\t", stringsAsFactors = F)
  silent.genes <- silent.genes$V1
  
  
  # 5. retrieve subsets of variable genes for mark 
  
  # 5.1. variable also for expression
  x <- intersect(variable.genes.expression, mark.sig.genes)
  
  # 5.2. silent for expression
  y <- intersect(silent.genes, mark.sig.genes)
  
  # 5.3. flat for expression
  z <- setdiff(mark.sig.genes, c(x, y))
  
  
  # 6. compute fold-change
  mark.m <- read.table(paste0(mark, "/QN.merged/", mark, ".matrix.after.QN.merged.tsv"), h=T, sep="\t")
  mark.m$FC <- apply(mark.m, 1, function(x){x <- max(x) - min(x)})
  
  
  # 7. prepare df for plot
  x.m <- mark.m[rownames(mark.m) %in% x, ]
  y.m <- mark.m[rownames(mark.m) %in% y, ]
  z.m <- mark.m[rownames(mark.m) %in% z, ]
  
  
  x.m$type <- paste0("expressed\nvariable (", nrow(x.m), ")")
  y.m$type <- paste0("silent (", nrow(y.m), ")")
  z.m$type <- paste0("expressed\nflat (", nrow(z.m), ")")
  
  df.plot <- rbind(x.m, y.m, z.m)
  df.plot$type <- factor(df.plot$type, levels = c(paste0("expressed\nvariable (", nrow(x.m), ")"),
                                                  paste0("expressed\nflat (", nrow(z.m), ")"), 
                                                  paste0("silent (", nrow(y.m), ")")))
  
  p <- ggplot(df.plot, aes(x=type, y=FC)) +
    geom_violin(aes(fill = type), trim = FALSE, color="white") +
    geom_boxplot(width = 0.2)+
    stat_compare_means( comparisons = list( c(paste0("expressed\nvariable (", nrow(x.m), ")"), 
                                              paste0("expressed\nflat (", nrow(z.m), ")")), 
                                            c(paste0("expressed\nflat (", nrow(z.m), ")"), 
                                              paste0("silent (", nrow(y.m), ")")), 
                                            c(paste0("expressed\nvariable (", nrow(x.m), ")"), 
                                              paste0("silent (", nrow(y.m), ")")) ) ) +
    labs(title = mark) +
    theme_bw() +  
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=15),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=15),
          strip.text.x = element_text(size=15),
          panel.border = element_rect(color="black"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          plot.title = element_text(size=15, hjust = .5),
          axis.line = element_line(colour = "black")) +
    expand_limits(y = 0) +
    scale_fill_manual(values = c("#525252", "#737373", "#bdbdbd")) +
    guides(fill=F)
  
  return(p)
  
}


#********
# BEGIN *
#********

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/")

marks <- c("H3K4me1", "H3K4me2", "H3K4me3", 
           "H3K9ac", "H3K27ac", "H4K20me1", 
           "H3K36me3", "H3K27me3", "H3K9me3")

lop <- list()
lop.peaks <- list()

for ( i in 1:9 ) {
  
  lop[[i]] <- my.function(mark = marks[i])
  lop.peaks[[i]] <- my.function(mark = marks[i], peaks=T)
  
}


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1z.pdf",
    width = 14, height = 14)
plot_grid(plotlist = lop, nrow=3, ncol=3, align="hv")
dev.off()


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1z.peaks.pdf",
    width = 14, height = 14)
plot_grid(plotlist = lop.peaks, nrow=3, ncol=3, align="hv")
dev.off()