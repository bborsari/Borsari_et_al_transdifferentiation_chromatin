.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")



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
library(stringr)



#********
# BEGIN *
#********


# 1. set working directory
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


# 2. the marks we're analyzing
marks <- c("H3K27ac", "H3K9ac", "H4K20me1", "H3K4me3", "H3K4me1",
           "H3K36me3", "H3K4me2", "H3K9me3", "H3K27me3")

# 3. read  ".6.groups.tsv" of each mark
# and store them in a list
lom <- list()
for ( i in 1:9 ) {
  
  tmp <- read.table(paste0(marks[i], "/QN.merged/", marks[i], ".6.groups.tsv"), 
                    h=T, sep="\t", stringsAsFactors = F)
  tmp$group <- gsub("peak_not_TSS", "no_peak", tmp$group) # merge peak_not_TSS and no_peak
  lom[[i]] <- tmp
  
}
rm(tmp)


# 4. check order of rownames
for (i in 2:9){
  stopifnot(identical(rownames(lom[[1]]), 
                      rownames(lom[[i]])))
}


# 5. prepare merged data.frame of groups across all marks
x <- data.frame(H3K27ac = lom[[1]]$group)

for (i in 2:9) {
  
  tmp <- data.frame(lom[[i]]$group)
  colnames(tmp) <- marks[i]
  x <- cbind(x, tmp)
  
}

rm(tmp)
rownames(x) <- rownames(lom[[1]])


# 6. cluster correspondence analysis
# performed in the cluster with R-3.5.2
# with the command below

# out.mca <- clusmca(x, nclus = 3, ndim = 3, method=c("clusCA","iFCB","MCAk"),
#         alphak = .5, nstart = 100, smartStart = NULL, gamma = TRUE,
#         seed = 1234)


load(".mca.RData")




# 5. dataframe with 3d coordinates
# of the attributes (combinations of patterns and marks)

## 5.1. extract attributes coordinates
mca_vars_df <- out.mca$attcoord
colnames(mca_vars_df) <- c("X", "Y", "Z")
## 5.2. add info of mark and group
mca_vars_df$mark <- str_split_fixed(rownames(mca_vars_df), "\\.", 2)[, 1]
mca_vars_df$group <- str_split_fixed(rownames(mca_vars_df), "\\.", 2)[, 2]
mca_vars_df$Z <- (mca_vars_df$Z)*(-1)
mca_vars_df$Y <- (mca_vars_df$Y)*(-1)


# 6. define shapes and color palette
shapes <- c("H3K4me1" = 0,
            "H3K4me2" = 7,
            "H3K27ac" = 1,
            "H3K9ac" = 10,
            "H3K4me3" = 2,
            "H3K36me3" = 5,
            "H4K20me1" = 9,
            "H3K9me3" = 4,
            "H3K27me3" = 8)

palette <- c("positively_correlated" = "#5B9F80",
             "no_peak" = "#000000",
             "negatively_correlated" = "#993404",
             "stable" = "#fec44f",
             "not_correlated" = "#ec7014")


# 7. 3d plot of the attributes 
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3g.marks.pdf", 
    width=5, height=4)
s3d <- scatterplot3d(mca_vars_df[mca_vars_df$group == "negatively_correlated", c(1, 3, 2)],
                     pch=shapes[mca_vars_df[mca_vars_df$group == "negatively_correlated", "mark"]],
                     color=palette[mca_vars_df[mca_vars_df$group == "negatively_correlated", "group"]],
                     angle = 40, grid = T, box = T,
                     cex.symbols = 1,
                     scale.y = 0.5,
                     xlim = c(-2, 14),
                     ylim = c(-5, 3),
                     zlim = c(-3, 2))
text(s3d$xyz.convert(mca_vars_df[mca_vars_df$group == "negatively_correlated", c(1, 3, 2)]), 
     labels = mca_vars_df[mca_vars_df$group == "negatively_correlated", "mark"])
dev.off()



# 8. dataframe with 3d coordinates
# of the objects (genes)
mca_genes_df <- out.mca$obscoord
colnames(mca_genes_df) <- c("X", "Y", "Z")
mca_genes_df$Z <- (mca_genes_df$Z)*(-1)
mca_genes_df$Y <- (mca_genes_df$Y)*(-1)
rownames(mca_genes_df) <- rownames(x)
mca_genes_df$cluster <- out.mca$cluster
mca_genes_df$cluster <- as.factor(mca_genes_df$cluster)
mca_genes_df$cluster <- paste0("cluster", mca_genes_df$cluster)


# 9. define color palette for clusters
palette2 <- c("cluster1" = "#F56E47", "cluster2" = "#97BF04", "cluster3" = "#772B59")


# 7. 3d plot of the objects (genes) 
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3g.genes.pdf", 
    width=5, height=4)
scatterplot3d(mca_genes_df[, c(1, 3, 2)], pch = 16, 
              color=palette2[mca_genes_df$cluster],
              angle = 40, grid = T, box = T,
              cex.symbols = 0.8,
              scale.y = 0.5,
              xlim = c(-2, 14),
              ylim = c(-5, 3),
              zlim = c(-3, 2))
dev.off()






# 8. retrieve clusters
# cl1 <- as.data.frame(rownames(mca_genes_df[mca_genes_df$cluster==1, ]))
# cl2 <- as.data.frame(rownames(mca_genes_df[mca_genes_df$cluster==2, ]))
# cl3 <- as.data.frame(rownames(mca_genes_df[mca_genes_df$cluster==3, ]))

# write.table(cl1, "~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.1.txt", row.names = F, col.names = F, quote=F)
# write.table(cl2, "~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.2.txt", row.names = F, col.names = F, quote=F)
# write.table(cl3, "~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.3.txt", row.names = F, col.names = F, quote=F)




