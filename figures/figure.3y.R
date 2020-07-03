.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")

palette <- c("positively_correlated" = "#5B9F80",
             "no_peak" = "#403734",
             "negatively_correlated" = "#993404",
             "stable" = "#fec44f",
             "not_correlated" = "#ec7014")



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

#********
# BEGIN *
#********

palette <- c("positively_correlated" = "#DEA450",
             "no_peak" = "#9B461F",
             "negatively_correlated" = "#a095a0",
             "stable" = "#5d6f62",
             "not_correlated" = "#42354C")


# 1. read marks 6 groups dataframes
marks <- c("H3K27ac", "H3K9ac", "H4K20me1", "H3K4me3", "H3K4me1",
           "H3K36me3", "H3K4me2", "H3K9me3", "H3K27me3")

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
cl1 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.1.txt")
cl1.groups <- all.marks.groups[rownames(all.marks.groups) %in% cl1$V1, ]
cl1.groups$gene_id <- rownames(cl1.groups)
cl1.groups[, 1:9] <- apply(cl1.groups[, 1:9], 2, as.character)

my.df <- data.frame(stringsAsFactors = F)
for ( i in 1:9) {
  
  for (j in i:9) {
    
    if (i!=j) {
      
      tmp <- melt(table(cl1.groups[, c(i,j)]))
      colnames(tmp) <- c("from", "to", "n")
      my.df <- rbind(my.df, tmp)
    }
    
  }
}


my.df$from_to <- paste(my.df$from, my.df$to, sep=";")
my.df.summary <- ddply(my.df,"from_to",numcolwise(sum))
my.df.summary$from <- str_split_fixed(my.df.summary$from_to, ";", 2)[, 1]
my.df.summary$to <- str_split_fixed(my.df.summary$from_to, ";", 2)[, 2]

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3y.cluster1.pdf", height=4.5, width=4.5)
chordDiagram(my.df.summary[, c("from", "to", "n")], order=names(palette), grid.col = palette)
dev.off()

# 4.2. - cluster 2
cl2 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.2.txt")
cl2.groups <- all.marks.groups[rownames(all.marks.groups) %in% cl2$V1, ]
cl2.groups$gene_id <- rownames(cl2.groups)
cl2.groups[, 1:9] <- apply(cl2.groups[, 1:9], 2, as.character)

my.df <- data.frame(stringsAsFactors = F)
for ( i in 1:9) {
  
  for (j in i:9) {
    
    if (i!=j) {
      
      tmp <- melt(table(cl2.groups[, c(i,j)]))
      colnames(tmp) <- c("from", "to", "n")
      my.df <- rbind(my.df, tmp)
    }
    
  }
}


my.df$from_to <- paste(my.df$from, my.df$to, sep=";")
my.df.summary <- ddply(my.df,"from_to",numcolwise(sum))
my.df.summary$from <- str_split_fixed(my.df.summary$from_to, ";", 2)[, 1]
my.df.summary$to <- str_split_fixed(my.df.summary$from_to, ";", 2)[, 2]

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3y.cluster2.pdf", height=4.5, width=4.5)
chordDiagram(my.df.summary[, c("from", "to", "n")], order=names(palette), grid.col = palette)
dev.off()


# 4.3. - cluster 3
cl3 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.3.txt")
cl3.groups <- all.marks.groups[rownames(all.marks.groups) %in% cl3$V1, ]
cl3.groups$gene_id <- rownames(cl3.groups)
cl3.groups[, 1:9] <- apply(cl3.groups[, 1:9], 2, as.character)

my.df <- data.frame(stringsAsFactors = F)
for ( i in 1:9) {
  
  for (j in i:9) {
    
    if (i!=j) {
      
      tmp <- melt(table(cl3.groups[, c(i,j)]))
      colnames(tmp) <- c("from", "to", "n")
      my.df <- rbind(my.df, tmp)
    }
    
  }
}


my.df$from_to <- paste(my.df$from, my.df$to, sep=";")
my.df.summary <- ddply(my.df,"from_to",numcolwise(sum))
my.df.summary$from <- str_split_fixed(my.df.summary$from_to, ";", 2)[, 1]
my.df.summary$to <- str_split_fixed(my.df.summary$from_to, ";", 2)[, 2]

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3y.cluster3.pdf", height=4.5, width=4.5)
chordDiagram(my.df.summary[, c("from", "to", "n")], order=names(palette), grid.col = palette)
dev.off()
