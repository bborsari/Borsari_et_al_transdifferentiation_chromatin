.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")



#************
# LIBRARIES *
#************


library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(circlize)
library(seriation)


#************
# FUNCTIONS *
#************


my.function1 <- function(mark) {
  
  # 1. read mark matrix of all genes after QN merged normalization
  mark.matrix <- read.table(paste0(mark, "/QN.merged/", mark, ".matrix.after.QN.merged.tsv"),
                            h=T, sep="\t")
  
  # 2. keep only flat genes (aka genes in expression matrix)
  mark.matrix <- mark.matrix[rownames(mark.matrix) %in% rownames(expression.matrix), ]
  
  # 3. check order of rows is the same as with expression
  mark.matrix <- mark.matrix[rownames(expression.matrix), ]
  stopifnot(identical(rownames(expression.matrix), rownames(mark.matrix)))
  
  # 4. center and scale mark matrix of flat genes 
  mark.matrix.scaled <- as.data.frame(t(scale(t(mark.matrix))))
  
  return(mark.matrix.scaled)
  
}




my.function2 <- function(k, mark) {
  
  # 1. read mark matrix of all genes after QN merged normalization
  mark.matrix.scaled <- marks.matrices[[k]]
    
  # 2. retrieve genes that are significantly variable for the mark
  mark.sig.genes <- read.table(paste0(mark, "/QN.merged/", mark, ".QN.merged.maSigPro.out.tsv"), h=T, sep="\t")
  mark.sig.genes$gene_id <- rownames(mark.sig.genes)
  mark.sig.genes <- mark.sig.genes %>% separate(gene_id, c("gene", "id"), "\\.")
  rownames(mark.sig.genes) <- mark.sig.genes$gene
  mark.sig.genes$gene <- NULL
  mark.sig.genes$id <- NULL
  mark.sig.genes <- rownames(mark.sig.genes)[rownames(mark.sig.genes) %in% rownames(expression.matrix)]
  
  # 3. retrieve genes that are flat for the mark
  mark.ns.genes <- setdiff(rownames(expression.matrix), mark.sig.genes)
  
  # 7. prepare separate expression and mark matrices for genes that are 
  # flat for expression only (sig)
  # flat for both expression and mark (ns)
  mark.matrix.scaled.sig <- mark.matrix.scaled[rownames(mark.matrix.scaled) %in% mark.sig.genes, ]

  # 8. retrieve row order for sig. genes with clustering
  ht.sig <- Heatmap(mark.matrix.scaled.sig,
                    col = palette2, 
                    name = mark, column_title = mark,
                    column_title_gp = gpar(fontsize = 60),
                    column_names_gp = gpar(fontsize = 60),
                    show_row_names = FALSE, width = unit(60, "mm"),
                    cluster_rows = T, cluster_columns = F,
                    clustering_distance_rows = "spearman",
                    clustering_method_rows = "complete",
                    show_column_names = F,
                    show_heatmap_legend = F,
                    heatmap_legend_param = list(title = "z-score"),
                    gap = unit(5, "mm"))
  
  new.order.sig <- rownames(mark.matrix.scaled.sig)[row_order(ht.sig)[[1]]]
  # mark.matrix.scaled.sig <- mark.matrix.scaled.sig[new.order.sig, ]
  # expression.matrix.scaled.sig <- expression.matrix.scaled.sig[new.order.sig, ]
  
  # new.order.ns <- rownames(mark.matrix.scaled.ns)
  # mark.matrix.scaled.ns <- mark.matrix.scaled.ns[new.order.ns, ]
  # expression.matrix.scaled.ns <- expression.matrix.scaled.ns[new.order.ns, ]
  
  
  # 10. reorder matrices for other marks
  
  for ( i in 1:9 ){
    
    tmp <- marks.matrices[[i]]
    tmp <- tmp[rownames(tmp) %in% new.order.sig, ]
    marks.matrices[[i]] <- tmp
    
  }
  
  
  return(marks.matrices)
  
  
  
  # # 10. plot the heatmap
  # mark_ht_list =
  #   
  #   Heatmap(marks.matrices[[1]],
  #           col = palette2, 
  #           name = "H3K27ac", column_title = "H3K27ac",
  #           column_title_gp = gpar(fontsize = 22),
  #           column_names_gp = gpar(fontsize = 22),
  #           show_row_names = FALSE, width = unit(40, "mm"),
  #           cluster_rows = F, cluster_columns = F,
  #           show_column_names = F,
  #           show_heatmap_legend = F,
  #           heatmap_legend_param = list(title = "z-score"),
  #           gap = unit(5, "mm")) +
  #   
  #   Heatmap(marks.matrices[[2]],
  #           col = palette2, 
  #           name = "H3K9ac", column_title = "H3K9ac",
  #           column_title_gp = gpar(fontsize = 22),
  #           column_names_gp = gpar(fontsize = 22),
  #           show_row_names = FALSE, width = unit(40, "mm"),
  #           cluster_rows = F, cluster_columns = F,
  #           show_column_names = F,
  #           show_heatmap_legend = F,
  #           heatmap_legend_param = list(title = "z-score"),
  #           gap = unit(5, "mm")) +
  #   
  #   Heatmap(marks.matrices[[3]],
  #           col = palette2, 
  #           name = "H4K20me1", column_title = "H4K20me1",
  #           column_title_gp = gpar(fontsize = 22),
  #           column_names_gp = gpar(fontsize = 22),
  #           show_row_names = FALSE, width = unit(40, "mm"),
  #           cluster_rows = F, cluster_columns = F,
  #           show_column_names = F,
  #           show_heatmap_legend = F,
  #           heatmap_legend_param = list(title = "z-score"),
  #           gap = unit(5, "mm")) +
  #   
  #   Heatmap(marks.matrices[[4]],
  #           col = palette2, 
  #           name = "H3K4me1", column_title = "H3K4me1",
  #           column_title_gp = gpar(fontsize = 22),
  #           column_names_gp = gpar(fontsize = 22),
  #           show_row_names = FALSE, width = unit(40, "mm"),
  #           cluster_rows = F, cluster_columns = F,
  #           show_column_names = F,
  #           show_heatmap_legend = F,
  #           heatmap_legend_param = list(title = "z-score"),
  #           gap = unit(5, "mm")) +
  #   
  #   Heatmap(marks.matrices[[5]],
  #           col = palette2, 
  #           name = "H3K4me3", column_title = "H3K4me3",
  #           column_title_gp = gpar(fontsize = 22),
  #           column_names_gp = gpar(fontsize = 22),
  #           show_row_names = FALSE, width = unit(40, "mm"),
  #           cluster_rows = F, cluster_columns = F,
  #           show_column_names = F,
  #           show_heatmap_legend = F,
  #           heatmap_legend_param = list(title = "z-score"),
  #           gap = unit(5, "mm")) +
  #   
  #   Heatmap(marks.matrices[[6]],
  #           col = palette2, 
  #           name = "H3K4me2", column_title = "H3K4me2",
  #           column_title_gp = gpar(fontsize = 22),
  #           column_names_gp = gpar(fontsize = 22),
  #           show_row_names = FALSE, width = unit(40, "mm"),
  #           cluster_rows = F, cluster_columns = F,
  #           show_column_names = F,
  #           show_heatmap_legend = F,
  #           heatmap_legend_param = list(title = "z-score"),
  #           gap = unit(5, "mm")) +
  #   
  #   Heatmap(marks.matrices[[7]],
  #           col = palette2, 
  #           name = "H3K9me3", column_title = "H3K9me3",
  #           column_title_gp = gpar(fontsize = 22),
  #           column_names_gp = gpar(fontsize = 22),
  #           show_row_names = FALSE, width = unit(40, "mm"),
  #           cluster_rows = F, cluster_columns = F,
  #           show_column_names = F,
  #           show_heatmap_legend = F,
  #           heatmap_legend_param = list(title = "z-score"),
  #           gap = unit(5, "mm")) +
  #   
  #   Heatmap(marks.matrices[[8]],
  #           col = palette2, 
  #           name = "H3K36me3", column_title = "H3K36me3",
  #           column_title_gp = gpar(fontsize = 22),
  #           column_names_gp = gpar(fontsize = 22),
  #           show_row_names = FALSE, width = unit(40, "mm"),
  #           cluster_rows = F, cluster_columns = F,
  #           show_column_names = F,
  #           show_heatmap_legend = F,
  #           heatmap_legend_param = list(title = "z-score"),
  #           gap = unit(5, "mm")) +
  #   
  #   Heatmap(marks.matrices[[9]],
  #           col = palette2, 
  #           name = "H3K27me3", column_title = "H3K27me3",
  #           column_title_gp = gpar(fontsize = 22),
  #           column_names_gp = gpar(fontsize = 22),
  #           show_row_names = FALSE, width = unit(40, "mm"),
  #           cluster_rows = F, cluster_columns = F,
  #           show_column_names = F,
  #           show_heatmap_legend = F,
  #           heatmap_legend_param = list(title = "z-score"),
  #           gap = unit(5, "mm"))
  # 
  # 
  # return(mark_ht_list)
  
}


#********
# BEGIN *
#********

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/")
palette2 = rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061'))


# 1. read expression matrix
expression.matrix <- read.table("expression/QN.merged/selected.genes.rep.2.3.after.QN.merged.tsv",
                                h=T, sep="\t")

# 2. retrieve list of silent genes
silent.genes <- read.table("expression/silent.genes.txt", h=F, sep="\t", stringsAsFactors = F)
silent.genes <- silent.genes$V1

# 3. retrieve list of variable genes
variable.genes <- read.table("expression/QN.merged/expression.matrix.tsv", h=T, sep="\t")
variable.genes <- rownames(variable.genes)

# 4. keep only flat genes 
expression.matrix <- expression.matrix[!(rownames(expression.matrix) %in%
                                           c(silent.genes, variable.genes)), ]

# 5. center and scale expression matrix of flat genes
expression.matrix.scaled <- as.data.frame(t(scale(t(expression.matrix))))


# 6. get matrices for all marks 

marks <- c("H3K27ac", "H3K9ac", "H4K20me1", 
           "H3K4me1", "H3K4me3", "H3K4me2", 
           "H3K9me3", "H3K36me3", "H3K27me3")
marks.matrices <- list()

for ( i in 1:9 ) {
  
  marks.matrices[[i]] <- my.function1(mark = marks[i])
  
}


marks.matrices <- my.function2(mark="H3K27ac", k=1)


# 7. compute correlations
marks.correlations <- data.frame(gene_id = rownames(marks.matrices[[1]]))

for ( i in 1:8 ) {

  for ( j in (i+1):9) {

    tmp <- diag(cor(t(marks.matrices[[i]]),
                    t(marks.matrices[[j]])))

    names(tmp) <- rownames(marks.matrices[[i]])

    marks.correlations <- cbind(marks.correlations, tmp)
    colnames(marks.correlations)[ncol(marks.correlations)] <- paste(marks[i],
                                                                    marks[j],
                                                                    sep = "_")

  }

}


marks.correlations$gene_id <- NULL
marks.correlations$median <- apply(marks.correlations, 1, median)
corr.order <- rownames(marks.correlations[order(marks.correlations$median,
                                               decreasing = T), ])
corr.order <- data.frame(gene_id = corr.order,
                         rank_corr = 1:1459)
cluster.order <- data.frame(gene_id = rownames(marks.matrices[[1]]),
                            rank_cluster = 1:1459)
corr.cluster.order <- merge(corr.order, cluster.order, by = "gene_id")
pheatmap(marks.matrices[[1]], cluster_cols = F, 
         clustering_distance_rows = "correlation",
         show_rownames = F,
         color = palette2)


test <- as.matrix(marks.matrices[[1]])

for ( i in 2:9 ) {
  
  test <- cbind(test, as.matrix(marks.matrices[[i]]))
  
}

test <- test[complete.cases(test), ]


o <- seriate(test)


Heatmap(test,
        row_order = get_order(o, 1),
        cluster_rows = F, cluster_columns = F)


# test <- marks.matrices[[1]]
# sig.genes <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/genes/flat.genes.expression/H3K27ac.txt",
#                         h=F, sep="\t", stringsAsFactors = F)
# sig.genes <- sig.genes$V1
# test <- test[rownames(test) %in% sig.genes, ]
# 
# test2 <- matrix(NA, nrow = 1459, ncol = 10)
# for ( i in 3:12 ) {
#   test2[, i-2] <- apply(test[, 1:i], 1, mean)
# }
# 
# test2 <- as.data.frame(test2)
# rownames(test2) <- rownames(test)
# 
# test2 <- test2[do.call(order, as.list(test2)),] 
# 
# Heatmap(test,
#         col = palette2,
#         cluster_rows = T, cluster_columns = F,
#         show_row_names = F, show_column_names = F,
#         clustering_distance_rows = "spearman")

