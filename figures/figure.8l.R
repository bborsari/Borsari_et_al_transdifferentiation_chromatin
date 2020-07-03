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
library(ggExtra)
library(ggalluvial)


#************
# FUNCTIONS *
#************

my.function <- function(x, y) {
  
  # 1. perform cluster correspondence analysis
  out.mca <- clusmca(all.marks.groups, nclus = x, ndim = y, method=c("clusCA","iFCB","MCAk"),
          alphak = .5, nstart = 100, smartStart = NULL, gamma = TRUE,
          seed = 1234)
  
  # 2. prepare df
  mca_genes_df  <- as.data.frame(out.mca$obscoord)
  rownames(mca_genes_df) <- rownames(all.marks.groups)

  set.seed(1234)
  # Empty vector to store the results
  res <- c()# Run K-means for different values of k and store J (WSS)
  for (k in 1:10){
    res <- rbind(res, cbind(k, J = kmeans(mca_genes_df, centers = k, nstart = 100)$tot.withinss))
  }
  p <- ggplot(as.data.frame(res), aes(k, J)) + 
    geom_path() + 
    geom_point() +
    scale_x_continuous(breaks = 1:10) +
    theme_bw()
  
  return(p)
  
}



#********
# BEGIN *
#********

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

# 5. best combination of number of clusters and dimensions
bestMethod <- tuneclus(data = all.marks.groups, 
                       nclusrange = 3:10, ndimrange = 2:9,
                       method = "clusCA",
                       nstart = 100, seed = 1234)

print(bestMethod)



# # 4. cluster correspondence analysis + elbow plot
# lop <- list()
# 
# lop[[1]] <- my.function(x = 3, y = 2)
# lop[[2]] <- my.function(x = 4, y = 2)
# lop[[3]] <- my.function(x = 5, y = 2)
# lop[[4]] <- my.function(x = 3, y = 3)
# lop[[5]] <- my.function(x = 4, y = 3)
# lop[[6]] <- my.function(x = 5, y = 3)
# lop[[7]] <- my.function(x = 4, y = 4)
# lop[[8]] <- my.function(x = 5, y = 4)
# 
# 
# pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8l.pdf",
#     width = 12, height=12)
# plot_grid(plotlist = lop, nrow = 3, ncol = 3,
#           align="hv")
# dev.off()
# 
# 
# 
# 
