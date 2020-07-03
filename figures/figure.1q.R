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

marks.matrices[[1]] <- expression.matrix.scaled

for ( i in 1:9 ) {
  
  marks.matrices[[i+1]] <- my.function1(mark = marks[i])
  
}


# 7. get input matrix for seriation 
seriate.m <- as.matrix(marks.matrices[[1]])

for ( i in 2:10 ) {
  
  seriate.m <- cbind(seriate.m, as.matrix(marks.matrices[[i]]))
  
}
seriate.m <- seriate.m[complete.cases(seriate.m), ]


# 8. heatmap
o <- seriate(seriate.m, method = "PCA_angle")
HT <-  Heatmap(seriate.m[, 1:12],
          row_order = get_order(o, 1),
          col = palette2, 
          name = "expression", column_title = "expression",
          column_title_gp = gpar(fontsize = 25),
          column_names_gp = gpar(fontsize = 25),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = F,
          heatmap_legend_param = list(title = "z-score"),
          gap = unit(5, "mm")) +
  
  Heatmap(seriate.m[, 13:24],
          row_order = get_order(o, 1),
          col = palette2, 
          name = "H3K27ac", column_title = "H3K27ac",
          column_title_gp = gpar(fontsize = 25),
          column_names_gp = gpar(fontsize = 25),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = F,
          heatmap_legend_param = list(title = "z-score"),
          gap = unit(5, "mm")) +
  
  Heatmap(seriate.m[, 25:36],
          row_order = get_order(o, 1),
          col = palette2, 
          name = "H3K9ac", column_title = "H3K9ac",
          column_title_gp = gpar(fontsize = 25),
          column_names_gp = gpar(fontsize = 25),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = F,
          heatmap_legend_param = list(title = "z-score"),
          gap = unit(5, "mm")) +
  
  Heatmap(seriate.m[, 37:48],
          row_order = get_order(o, 1),
          col = palette2, 
          name = "H4K20me1", column_title = "H4K20me1",
          column_title_gp = gpar(fontsize = 25),
          column_names_gp = gpar(fontsize = 25),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = F,
          heatmap_legend_param = list(title = "z-score"),
          gap = unit(5, "mm")) +
  
  Heatmap(seriate.m[, 49:60],
          row_order = get_order(o, 1),
          col = palette2, 
          name = "H3K4me1", column_title = "H3K4me1",
          column_title_gp = gpar(fontsize = 25),
          column_names_gp = gpar(fontsize = 25),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = F,
          heatmap_legend_param = list(title = "z-score"),
          gap = unit(5, "mm")) +
  
  Heatmap(seriate.m[, 61:72],
          row_order = get_order(o, 1),
          col = palette2, 
          name = "H3K4me3", column_title = "H3K4me3",
          column_title_gp = gpar(fontsize = 25),
          column_names_gp = gpar(fontsize = 25),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = F,
          heatmap_legend_param = list(title = "z-score"),
          gap = unit(5, "mm")) +
  
  Heatmap(seriate.m[, 73:84],
          row_order = get_order(o, 1),
          col = palette2, 
          name = "H3K4me2", column_title = "H3K4me2",
          column_title_gp = gpar(fontsize = 25),
          column_names_gp = gpar(fontsize = 25),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = F,
          heatmap_legend_param = list(title = "z-score"),
          gap = unit(5, "mm")) +
  
  Heatmap(seriate.m[, 85:96],
          row_order = get_order(o, 1),
          col = palette2, 
          name = "H3K9me3", column_title = "H3K9me3",
          column_title_gp = gpar(fontsize = 25),
          column_names_gp = gpar(fontsize = 25),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = F,
          heatmap_legend_param = list(title = "z-score"),
          gap = unit(5, "mm")) +
  
  Heatmap(seriate.m[, 97:108],
          row_order = get_order(o, 1),
          col = palette2, 
          name = "H3K36me3", column_title = "H3K36me3",
          column_title_gp = gpar(fontsize = 25),
          column_names_gp = gpar(fontsize = 25),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = F,
          heatmap_legend_param = list(title = "z-score"),
          gap = unit(5, "mm")) +
  
  Heatmap(seriate.m[, 109:120],
          row_order = get_order(o, 1),
          col = palette2, 
          name = "H3K27me3", column_title = "H3K27me3",
          column_title_gp = gpar(fontsize = 25),
          column_names_gp = gpar(fontsize = 25),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = F,
          heatmap_legend_param = list(title = "z-score"),
          gap = unit(5, "mm"))


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1q.pdf",
    width = 25, height = 6)
print(HT)
dev.off()
