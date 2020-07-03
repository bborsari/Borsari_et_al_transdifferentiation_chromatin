.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(pheatmap)
library(MASS)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(ggrepel)
library(plotly)
library(clustrd)
library(ggalluvial)
library(stringr)
library(circlize)

#************
# FUNCTIONS *
#************

my.function <- function(cluster) {
  
  cl <- read.table(paste0("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.", cluster, ".txt"))
  cl.groups <- all.marks.groups[rownames(all.marks.groups) %in% cl$V1, ]
  cl.groups$gene_id <- rownames(cl.groups)
  cl.groups[, 1:9] <- apply(cl.groups[, 1:9], 2, as.character)
  
  my.df <- data.frame(stringsAsFactors = F)
  
  for ( i in 1:9) {
    
    for (j in i:9) {
      
      if (i!=j) {
        
        tmp <- melt(table(cl.groups[, c(i,j)]))
        colnames(tmp) <- c("from", "to", "n")
        tmp$mark1 <- colnames(cl.groups)[i]
        tmp$mark2 <- colnames(cl.groups)[j]
        my.df <- rbind(my.df, tmp)
        
      }
      
    }
  }
  
  my.combined <- unique(data.frame(marks = c(my.df$mark1, my.df$mark2),
                                  groups = c(as.character(my.df$from), as.character(my.df$to)), 
                                  stringsAsFactors = FALSE))
  
  my.combined$groups <- factor(my.combined$groups, 
                               levels = c("no_peak",
                                          "positively_correlated",
                                          "stable",
                                          "not_correlated",
                                          "negatively_correlated"))
  
  my.combined$marks <- factor(my.combined$marks, levels = marks)
  
  
  my.combined <- my.combined[order(my.combined$groups, my.combined$marks), ]
  my.order <- paste(my.combined$groups, my.combined$marks, sep = "|")
  
  my.grid.col <- structure(palette[my.combined$groups], names = my.order)
  
  my.df2 = data.frame(from = paste(my.df$from, my.df$mark1, sep = "|"),
                      to = paste(my.df$to, my.df$mark2, sep = "|"),
                      value = my.df$n, stringsAsFactors = FALSE)
  
  my.gap = rep(1, length(my.order))
  my.gap[which(!duplicated(my.combined$groups, fromLast = TRUE))] = 5
  
  circos.clear()
  
  circos.par(gap.degree = my.gap)
  chordDiagram(my.df2, order = my.order, annotationTrack = c("grid", "axis"),
               grid.col = my.grid.col, directional = TRUE,
               preAllocateTracks = list(
                 track.height = 0.04,
                 track.margin = c(0.05, 0)
               )
  )
  
  
  for( mark in unique(my.combined$marks) ) {
    l = my.combined$marks == mark
    sn = paste(my.combined$groups[l], my.combined$marks[l], sep = "|")
    highlight.sector(sn, track.index = 1, col = palette2[mark], 
                     niceFacing = TRUE)
  }
  
  
  legend("bottomleft", pch = 15, col = palette, 
         legend = names(palette), cex = 0.6)
  legend("bottomright", pch = 15, col = palette2, 
         legend = names(palette2), cex = 0.6)
  
  circos.clear()
  
}



#********
# BEGIN *
#********

marks <- c("H3K27ac", "H3K9ac", "H4K20me1", "H3K4me3", "H3K4me1",
           "H3K36me3", "H3K4me2", "H3K9me3", "H3K27me3")


palette <- c("no_peak" = "#403734",
             "positively_correlated" = "#5B9F80",
             "stable" = "#fec44f",
             "not_correlated" = "#ec7014",
             "negatively_correlated" = "#993404")


palette2 <- c("H3K9ac" = "#af4d85",
              "H3K27ac" = "#630039",
              "H3K4me3" = "#d199b9",
              "H3K27me3" = "#1d2976",
              "H3K9me3" = "#a7add4",
              "H3K36me3" = "#7fbc41",
              "H4K20me1" = "#4c7027",
              "H3K4me1" = "#e5ab00",
              "H3K4me2" = "#a67c00")



# 1. read marks 6 groups dataframes

all.marks.groups.list <- list()

for ( i in 1:9 ) {
  
  tmp <- read.table(paste0(marks[i], "/QN.merged/", marks[i], ".6.groups.tsv"), h=T, sep="\t", stringsAsFactors = F)
  tmp$group <- gsub("peak_not_TSS", "no_peak", tmp$group)
  tmp$final_class <- gsub("regulation", "-\nregulated", tmp$final_class)
  all.marks.groups.list[[i]] <- tmp
  
}

rm(tmp)


# 2. check order of rownames
for (i in 2:9){
  
  stopifnot(identical(rownames(all.marks.groups.list[[1]]), 
                      rownames(all.marks.groups.list[[i]])))
}


# 3. prepare merged data.frame of groups across all marks
all.marks.groups <- data.frame(H3K27ac = all.marks.groups.list[[1]]$group)

for (i in 2:9) {
  
  tmp <- data.frame(all.marks.groups.list[[i]]$group)
  colnames(tmp) <- marks[i]
  all.marks.groups <- cbind(all.marks.groups, tmp)
  
}

rm(tmp)
rownames(all.marks.groups) <- rownames(all.marks.groups.list[[1]])


# 4. make circular plots

# 4.1. - cluster 1
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3y2.cluster1.pdf", 
    height=7, width=4.5)
print(my.function(cluster=1))
dev.off()


# 4.2. - cluster 2
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3y2.cluster2.pdf", 
    height=7, width=4.5)
print(my.function(cluster=2))
dev.off()


# 4.3. - cluster 3
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3y2.cluster3.pdf", 
    height=7, width=4.5)
print(my.function(cluster=3))
dev.off()

