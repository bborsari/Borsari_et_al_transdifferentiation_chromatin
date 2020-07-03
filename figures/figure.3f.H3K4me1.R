#************
# LIBRARIES *
#************

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)
library(circlize)


#************
# FUNCTIONS *
#************



col_fun = colorRamp2(c(-10, -1.74, -1, 2), # good
                     c("#c0dbeb", "#d9e9f3", "#f4a582", "#ca0020"))


heatmap.TSS.pc <- function(mark, class, my.title, group) {
  
  setwd(paste0("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/",
               mark, 
               "/QN.merged/bwtool.matrix/mark/"))
  
  m.H000X1 <- read.table(paste0("H000", mark, "X1.major.", group, ".", class, ".bwtool.matrix.tsv"), h=T, sep="\t", quote=NULL)
  
  m.H000X2 <- read.table(paste0("H000", mark, "X2.major.", group, ".", class, ".bwtool.matrix.tsv"), h=T, sep="\t", quote=NULL)
  
  m.H003X1 <- read.table(paste0("H003", mark, "X1.major.", group, ".", class, ".bwtool.matrix.tsv"), h=T, sep="\t", quote=NULL)
  
  m.H003X2 <- read.table(paste0("H003", mark, "X2.major.", group, ".", class, ".bwtool.matrix.tsv"), h=T, sep="\t", quote=NULL)
  
  m.H120X1 <- read.table(paste0("H120", mark, "X1.major.", group, ".", class, ".bwtool.matrix.tsv"), h=T, sep="\t", quote=NULL)
  
  m.H120X2 <- read.table(paste0("H120", mark, "X2.major.", group, ".", class, ".bwtool.matrix.tsv"), h=T, sep="\t", quote=NULL)
  
  m.H168X1 <- read.table(paste0("H168", mark, "X1.major.", group, ".", class, ".bwtool.matrix.tsv"), h=T, sep="\t", quote=NULL)
  
  m.H168X2 <- read.table(paste0("H168", mark, "X2.major.", group, ".", class, ".bwtool.matrix.tsv"), h=T, sep="\t", quote=NULL)
  
  stopifnot(identical(rownames(m.H000X1), rownames(m.H000X2)))
  stopifnot(identical(rownames(m.H000X1), rownames(m.H003X1)))
  stopifnot(identical(rownames(m.H000X1), rownames(m.H003X2)))
  stopifnot(identical(rownames(m.H000X1), rownames(m.H120X1)))
  stopifnot(identical(rownames(m.H000X1), rownames(m.H120X2)))
  stopifnot(identical(rownames(m.H000X1), rownames(m.H168X1)))
  stopifnot(identical(rownames(m.H168X1), rownames(m.H168X2)))
  
  
  m.H000X1.H000X2 <- (m.H000X1 + m.H000X2) / 2
  m.H003X1.H003X2 <- (m.H003X1 + m.H003X2) / 2
  m.beginning <- (m.H000X1.H000X2 + m.H003X1.H003X2) / 2
  
  m.H120X1.H120X2 <- (m.H120X1 + m.H120X2) / 2
  m.H168X1.H168X2 <- (m.H168X1 + m.H168X2) / 2
  m.end <- (m.H120X1.H120X2 + m.H168X1.H168X2) / 2
  
  
  mat <- cbind(m.beginning,
               m.end)
  
  
  mat.mean <- data.frame(my.mean1= apply(mat, 1, mean))
  
  mat.mean <- mat.mean[order(mat.mean$my.mean1,
                             decreasing = T), , drop=F]
  
  mat <- mat[rownames(mat.mean), ]
  
  
  p =  Heatmap(log2(mat + 0.001),
               cluster_rows = F,
               cluster_columns = F,
               show_heatmap_legend = F,
               show_row_names = FALSE, 
               width = unit(100, "mm"),
               show_column_names = F,
               col = col_fun,
               column_split = rep(c("0h - 3h", "120h - 168h"), each  = 4000),
               border = T)
  
  p <- draw(p, column_title = my.title)
  
  return(p)
  
}

heatmap.TSS.stable <- function(mark, my.title, group) {
  
  setwd(paste0("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/",
               mark, 
               "/QN.merged/bwtool.matrix/mark/"))
  
  
  # H000X1
  m.H000X1.upreg <- read.table(paste0("H000", mark, "X1.major.", group, ".upregulation.bwtool.matrix.tsv"), 
                               h=T, sep="\t", quote=NULL)
  m.H000X1.downreg <- read.table(paste0("H000", mark, "X1.major.", group, ".downregulation.bwtool.matrix.tsv"), 
                                 h=T, sep="\t", quote=NULL)
  m.H000X1 <- rbind(m.H000X1.upreg, m.H000X1.downreg)
  
  
  # H000X2
  m.H000X2.upreg <- read.table(paste0("H000", mark, "X2.major.", group, ".upregulation.bwtool.matrix.tsv"), 
                               h=T, sep="\t", quote=NULL)
  m.H000X2.downreg <- read.table(paste0("H000", mark, "X2.major.", group, ".downregulation.bwtool.matrix.tsv"), 
                                 h=T, sep="\t", quote=NULL)
  m.H000X2 <- rbind(m.H000X2.upreg, m.H000X2.downreg)
  
  
  # H003X1
  m.H003X1.upreg <- read.table(paste0("H003", mark, "X1.major.", group, ".upregulation.bwtool.matrix.tsv"), 
                               h=T, sep="\t", quote=NULL)
  m.H003X1.downreg <- read.table(paste0("H003", mark, "X1.major.", group, ".downregulation.bwtool.matrix.tsv"), 
                                 h=T, sep="\t", quote=NULL)
  m.H003X1 <- rbind(m.H003X1.upreg, m.H003X1.downreg)
  
  
  # H003X2
  m.H003X2.upreg <- read.table(paste0("H003", mark, "X2.major.", group, ".upregulation.bwtool.matrix.tsv"), 
                               h=T, sep="\t", quote=NULL)
  m.H003X2.downreg <- read.table(paste0("H003", mark, "X2.major.", group, ".downregulation.bwtool.matrix.tsv"), 
                                 h=T, sep="\t", quote=NULL)
  m.H003X2 <- rbind(m.H003X2.upreg, m.H003X2.downreg)
  
  
  # H120X1
  m.H120X1.upreg <- read.table(paste0("H120", mark, "X1.major.", group, ".upregulation.bwtool.matrix.tsv"), 
                               h=T, sep="\t", quote=NULL)
  m.H120X1.downreg <- read.table(paste0("H120", mark, "X1.major.", group, ".downregulation.bwtool.matrix.tsv"), 
                                 h=T, sep="\t", quote=NULL)
  m.H120X1 <- rbind(m.H120X1.upreg, m.H120X1.downreg)
  
  
  # H120X2
  m.H120X2.upreg <- read.table(paste0("H120", mark, "X2.major.", group, ".upregulation.bwtool.matrix.tsv"), 
                               h=T, sep="\t", quote=NULL)
  m.H120X2.downreg <- read.table(paste0("H120", mark, "X2.major.", group, ".downregulation.bwtool.matrix.tsv"), 
                                 h=T, sep="\t", quote=NULL)
  m.H120X2 <- rbind(m.H120X2.upreg, m.H120X2.downreg)
  
  
  m.H168X1.upreg <- read.table(paste0("H168", mark, "X1.major.", group, ".upregulation.bwtool.matrix.tsv"), 
                               h=T, sep="\t", quote=NULL)
  m.H168X1.downreg <- read.table(paste0("H168", mark, "X1.major.", group, ".downregulation.bwtool.matrix.tsv"), 
                                 h=T, sep="\t", quote=NULL)
  m.H168X1 <- rbind(m.H168X1.upreg, m.H168X1.downreg)
  
  
  # H168X2
  m.H168X2.upreg <- read.table(paste0("H168", mark, "X2.major.", group, ".upregulation.bwtool.matrix.tsv"), 
                               h=T, sep="\t", quote=NULL)
  m.H168X2.downreg <- read.table(paste0("H168", mark, "X2.major.", group, ".downregulation.bwtool.matrix.tsv"), 
                                 h=T, sep="\t", quote=NULL)
  m.H168X2 <- rbind(m.H168X2.upreg, m.H168X2.downreg)
  
  stopifnot(identical(rownames(m.H000X1), rownames(m.H000X2)))
  stopifnot(identical(rownames(m.H000X1), rownames(m.H003X1)))
  stopifnot(identical(rownames(m.H000X1), rownames(m.H003X2)))
  stopifnot(identical(rownames(m.H000X1), rownames(m.H120X1)))
  stopifnot(identical(rownames(m.H000X1), rownames(m.H120X2)))
  stopifnot(identical(rownames(m.H000X1), rownames(m.H168X1)))
  stopifnot(identical(rownames(m.H168X1), rownames(m.H168X2)))
  
  
  m.H000X1.H000X2 <- (m.H000X1 + m.H000X2) / 2
  m.H003X1.H003X2 <- (m.H003X1 + m.H003X2) / 2
  m.beginning <- (m.H000X1.H000X2 + m.H003X1.H003X2) / 2
  
  m.H120X1.H120X2 <- (m.H120X1 + m.H120X2) / 2
  m.H168X1.H168X2 <- (m.H168X1 + m.H168X2) / 2
  m.end <- (m.H120X1.H120X2 + m.H168X1.H168X2) / 2
  
  mat <- cbind(m.beginning,
               m.end)
  
  mat.mean <- data.frame(my.mean1= apply(mat, 1, mean))
  
  mat.mean <- mat.mean[order(mat.mean$my.mean1,
                             decreasing = T), , drop=F]
  
  mat <- mat[rownames(mat.mean), ]
  
  
  p =  Heatmap(log2(mat + 0.001),
               cluster_rows = F,
               cluster_columns = F,
               show_heatmap_legend = F,
               show_row_names = FALSE, 
               width = unit(100, "mm"),
               show_column_names = F,
               col = col_fun,
               column_split = rep(c("0h - 3h", "120h - 168h"), each  = 4000),
               border = T)
  
  p <- draw(p, column_title = my.title)
  
  return(p)
  
}



#********
# BEGIN *
#********


# 1. upreg & positively_correlated
p.H3K4me1.upregulation.pc <- heatmap.TSS.pc(mark = "H3K4me1",
                                            class = "upregulation",
                                            group = "positively_correlated",
                                            my.title = "H3K4me1 - upreg. - pc")

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3f.H3K4me1.upregulation.positively_correlated.pdf",
    width = 5,
    height = 5)
p.H3K4me1.upregulation.pc
dev.off()


# 2. downreg & positively_correlated
p.H3K4me1.downregulation.pc <- heatmap.TSS.pc(mark = "H3K4me1",
                                              class = "downregulation",
                                              group = "positively_correlated",
                                              my.title = "H3K4me1 - downreg. - pc")

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3f.H3K4me1.downregulation.positively_correlated.pdf",
    width = 5,
    height = 5)
p.H3K4me1.downregulation.pc
dev.off()


# 3. upreg + downreg & stable
p.H3K4me1.stable <- heatmap.TSS.stable(mark = "H3K4me1",
                                       group = "stable",
                                       my.title = "H3K4me1 - stable")


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3f.H3K4me1.stable.pdf",
    width = 5,
    height = 5)
p.H3K4me1.stable
dev.off()
