.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/expression/QN.merged")


#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(pheatmap)


#********
# BEGIN *
#********

# 1. read dataframes
expression.matrix <- read.table("selected.genes.rep.2.3.after.QN.merged.tsv", h=T, sep="\t")
metadata <- read.table("metadata.tsv", h=T, sep="\t")
metadata$gene_id <- rownames(metadata)


# 2. perform clustering of high_expression, moderate_expression, low_expression, variable
# according to expression

# 2.1. retrieve expression matrix of genes that failed first step of classification
expression.matrix.subset <- expression.matrix[rownames(expression.matrix) %in%
                                                rownames(metadata[metadata$class %in% c("high_expression",
                                                                                        "moderate_expression",
                                                                                        "variable"), ]), ]
# 2.2. heatmap for exploratory analysis
pheatmap(t(scale(t(expression.matrix.subset))),
         show_rownames = F,
         cluster_cols = F,
         clustering_distance_rows = "euclidean",
         clustering_method = "complete")

# 2.3. perform hierarchical clustering
X <- dist(t(scale(t(expression.matrix.subset)))) # distance = euclidean
hca <- hclust(X) # method = complete
# new.order <- hca$order
# my.genes <- rownames(metadata[metadata$class %in% c("high_expression", 
#                                                     "moderate_expression",
#                                                     "variable"), ])[new.order]


# 2.4. cut the tree in 3 clusters (up-, down-regulated and bending)
further_classification <- data.frame(cutree(hca, k=3))
colnames(further_classification) <- "hc"
further_classification$hc <- paste0("cluster_", further_classification$hc)
further_classification$gene_id <- rownames(further_classification)
class <- c()
for (r in rownames(further_classification)) {
  
  if (further_classification[r, "hc"] == "cluster_1") {
    
    tmp_class <- "downregulation"
    
  } else if (further_classification[r, "hc"] == "cluster_2") {
    
    tmp_class <- "upregulation"
    
  } else {
    
    tmp_class <- "bending"
    
  }
  
  class <- c(class, tmp_class)
}

# 2.5. re-plot results with cluster annotation
pheatmap(t(scale(t(expression.matrix.subset))),
         show_rownames = F,
         cluster_cols = F,
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         annotation_row = further_classification[, 1, drop=F])


# 2.6. add new classification to the metadata file
further_classification$class2 <- class
metadata <- merge( x = metadata,
                   y = further_classification,
                   by = "gene_id",
                   all = T)
rownames(metadata) <- metadata$gene_id

# 2.7. save new metadata file
write.table(metadata[, 2:ncol(metadata)], "metadata.class2.tsv", sep="\t", row.names = T, col.names = T, quote=F)

