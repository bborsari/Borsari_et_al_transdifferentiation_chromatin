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


heatmap.TSS <- function(mark, class, my.title) {
  
  groups <- c("positively_correlated", "stable", "no_peak")
  
  setwd(paste0("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/",
               mark, 
               "/QN.merged/bwtool.matrix/mark/"))
  
  this.mat <- data.frame()
  
  for (group in groups) {
    
    
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
    
    if (group == "positively_correlated") {
      
      aa = nrow(mat)
      
    } else if (group == "stable") {
      
      bb = nrow(mat)
      
    } else if (group == "negatively_correlated") {
      
      cc = nrow(mat)
      
    } else if (group == "no_peak") {
      
      dd = nrow(mat)
      
    } else {
      
      ee = nrow(mat)
      
    }
    
    
    mat.mean <- data.frame(my.mean1= apply(mat, 1, mean))
    
    mat.mean <- mat.mean[order(mat.mean$my.mean1,
                               decreasing = T), , drop=F]
    
    this.mat <- rbind(this.mat, mat[rownames(mat.mean), ])
    
  }
  
  # c("#ffffff", '#f0f0f0', '#f0f0f0', '#d9d9d9', 
  #   '#bdbdbd','#969696','#737373',
  #   '#525252', '#252525','#000000')
  

  # col_fun = colorRamp2(c(-10, 0, 0.59, 
  #                        1, 1.32, 1.59, 
  #                        1.81, 2, 2.17, 2.5), 
  #                      c("#053061", "#2166ac", "#4393c3", 
  #                        "#92c5de", "#d1e5f0", "#fddbc7", 
  #                        "#f4a582", "#d6604d","#b2182b", "#67001f"))
  
  # col_fun = colorRamp2(c(-10, -1, -0.7, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.35, 1.5, 1.6, 1.75, 2, 2.25, 2.4),
  #                      c("#081d58", "#253494", "#225ea8", "#1d91c0", "#41b6c4", "#7fcdbb", "#c7e9b4", "#edf8b1", "#ffffd9",
  #                        '#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026'))
  # 
  col_fun = colorRamp2(c(-10, 1, 1.8, 2.5), # good
                       c("#0571b0", "#92c5de", "#f4a582", "#ca0020"))
  

  
  p =  Heatmap(log2(this.mat + 0.001),
               cluster_rows = F,
               cluster_columns = F,
               show_heatmap_legend = T,
               show_row_names = FALSE, 
               width = unit(100, "mm"),
               show_column_names = F,
               col = col_fun,
               row_split = c(rep("1 - positively correlated", aa),
                             rep("2 - stable", bb),
                             rep("3 - no_peak", dd)),
                             # rep("3 - negatively_correlated", cc),
                             # rep("5 - not_correlated", ee)),
               column_split = rep(c("000h", "168h"), each  = 4000),
               border = T)
  
  p <- draw(p, column_title = my.title)

  return(p)

}

p.H3K4me3.upregulation <- heatmap.TSS(mark = "H3K4me3", class = "upregulation",
                                      my.title = "H3K4me3 - upreg.")
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3f2.H3K4me3.upregulation.pdf", width = 8)
p.H3K4me3.upregulation
dev.off()


# p.H3K4me3.downregulation <- heatmap.TSS(mark = "H3K4me3", class = "downregulation", j="TSS",
#                                         # palette =  c("#fff7f3", "#fff7f3", "#fff7f3", "#fff7f3", "#fff7f3",
#                                         #              "#fde0dd", "#fcc5c0", "#dd3497",
#                                         #              "#ae017e",
#                                         #              "#7a0177",
#                                         #              "#49006a"),
#                                         palette = c("#ffffff", "#ffffff", "#ffffff", "#ffffff", "#ffffff",
#                                                     '#f0f0f0', '#d9d9d9', '#737373', 
#                                                     '#525252', '#252525','#000000'),
#                                         my.title = "H3K4me3 - downreg.")
# 
# 
# pdf("~/public_html/paper_ERC/fig.3f.H3K4me3.downregulation.pdf", width = 5)
# p.H3K4me3.downregulation
# dev.off()
# 
# 
# p.H3K9ac.upregulation <- heatmap.TSS(mark = "H3K9ac", class = "upregulation", j = "TSS",
#                                      # palette = colorRampPalette(c("#fff7f3", "#fff7f3", "#fff7f3",
#                                      #                              "#fde0dd", "#fcc5c0", "#dd3497",
#                                      #                              "#ae017e",
#                                      #                              "#7a0177",
#                                      #                              "#49006a"))(100),
#                                      palette = colorRampPalette(c('#ffffff', '#ffffff', '#ffffff',
#                                                                   '#f0f0f0', '#d9d9d9', '#737373',
#                                                                   '#525252','#252525','#000000'))(100),
#                                      my.title = "H3K9ac - upreg.")
# pdf("~/public_html/paper_ERC/fig.3f.H3K9ac.upregulation.pdf", width = 5)
# p.H3K9ac.upregulation
# dev.off()
# 
# p.H3K9ac.downregulation <- heatmap.TSS(mark = "H3K9ac", class = "downregulation", j="TSS",
#                                        # palette = colorRampPalette(c("#fff7f3", "#fff7f3", "#fff7f3",
#                                        #                              "#fde0dd", "#fcc5c0", "#dd3497",
#                                        #                              "#ae017e",
#                                        #                              "#7a0177",
#                                        #                              "#49006a"))(100),
#                                        palette = colorRampPalette(c('#ffffff', '#ffffff', '#ffffff',
#                                                                     '#f0f0f0', '#d9d9d9', '#737373',
#                                                                     '#525252','#252525','#000000'))(100),
#                                        my.title = "H3K9ac - downreg.")
# 
# pdf("~/public_html/paper_ERC/fig.3f.H3K9ac.downregulation.pdf", width=5)
# p.H3K9ac.downregulation
# dev.off()
# 
# 
# p.H3K27ac.upregulation <- heatmap.TSS(mark = "H3K27ac", class = "upregulation", j="TSS",
#                                       # palette = c("#fff7f3", "#fff7f3",
#                                       #             "#fde0dd", "#fcc5c0", "#dd3497",
#                                       #             "#ae017e",
#                                       #             "#7a0177",
#                                       #             "#49006a"),
#                                       palette = c('#ffffff', '#ffffff',
#                                                   '#f0f0f0','#d9d9d9', '#737373',
#                                                   '#525252','#252525','#000000'),
#                                       my.title = "H3K27ac (TSS) - upreg.")
# 
# pdf("~/public_html/paper_ERC/fig.3f.H3K27ac_TSS.upregulation.pdf", width=5)
# p.H3K27ac.upregulation
# dev.off()
# 
# # ['#fff7f3','#fde0dd','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a']
# # ['#ffffff','#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525','#000000']
# 
# 
# p.H3K27ac.downregulation <- heatmap.TSS(mark = "H3K27ac", class = "downregulation", j="TSS",
#                                         # palette = c("#fff7f3",
#                                         #             "#fde0dd", "#fcc5c0", "#dd3497",
#                                         #             "#ae017e",
#                                         #             "#7a0177",
#                                         #             "#49006a"),
#                                         palette = c('#ffffff',
#                                                     '#f0f0f0','#d9d9d9', '#737373',
#                                                     '#525252','#252525','#000000'),
#                                         my.title = "H3K27ac (TSS) - downreg.")
# pdf("~/public_html/paper_ERC/fig.3f.H3K27ac_TSS.downregulation.pdf", width=5)
# p.H3K27ac.downregulation
# dev.off()
