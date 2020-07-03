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


my.function <- function(mark) {
  
  # 0. read list of genes with peaks in at least one time-point
  mark.genes.peaks <- read.table(paste0(mark, "/QN.merged/all.genes.intersecting.peaks.tsv"), 
                                 h=F, sep="\t", stringsAsFactors = F)
  mark.genes.peaks <- mark.genes.peaks$V1
  
  
  # 1. retrieve genes that are significantly variable for the mark
  # while silent for expression
  mark.sig.genes <- read.table(paste0(mark, "/QN.merged/", mark, ".QN.merged.maSigPro.out.tsv"), h=T, sep="\t")
  mark.sig.genes$gene_id <- rownames(mark.sig.genes)
  mark.sig.genes <- mark.sig.genes %>% separate(gene_id, c("gene", "id"), "\\.")
  rownames(mark.sig.genes) <- mark.sig.genes$gene
  mark.sig.genes$gene <- NULL
  mark.sig.genes$id <- NULL
  mark.sig.genes <- rownames(mark.sig.genes)[rownames(mark.sig.genes) %in% silent.genes]
  
  
  # 2. keep genes variable for the mark that also have a peak 
  mark.sig.genes <- intersect(mark.sig.genes, mark.genes.peaks)
  
  
  # 3. return genes that are not variable 
  mark.ns.genes <- setdiff(silent.genes, mark.sig.genes)
  return(mark.ns.genes)
  
}

my.function1 <- function(mark) {
  
  # 1. read mark matrix of all genes after QN merged normalization
  mark.matrix <- read.table(paste0(mark, "/QN.merged/", mark, ".matrix.after.QN.merged.tsv"),
                            h=T, sep="\t")
  
  # 2. keep only silent genes
  mark.matrix <- mark.matrix[rownames(mark.matrix) %in% silent.genes, ]
  
  # 3. check order of rows is the same as within silent.genes vector
  mark.matrix <- mark.matrix[silent.genes, ]
  stopifnot(identical(silent.genes, rownames(mark.matrix)))
  
  # 4. center and scale mark matrix of silent genes 
  mark.matrix.scaled <- as.data.frame(t(scale(t(mark.matrix))))
  
  return(mark.matrix.scaled)
  
}


#********
# BEGIN *
#********

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/")
palette2 = rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061'))


# 1. retrieve list of silent genes
silent.genes <- read.table("expression/silent.genes.txt", h=F, sep="\t", stringsAsFactors = F)
silent.genes <- silent.genes$V1


# 2. get matrices for all marks 
marks <- c("H3K4me1", "H3K4me2", "H3K4me3", 
           "H3K9ac", "H3K27ac", "H4K20me1", 
           "H3K36me3", "H3K27me3", "H3K9me3")
marks.matrices <- list()


for ( i in 1:9 ) {
  
  marks.matrices[[i]] <- my.function1(mark = marks[i])
  
}


# 3. get input matrix for seriation 
seriate.m <- as.matrix(marks.matrices[[1]])

for ( i in 2:9 ) {
  
  seriate.m <- cbind(seriate.m, as.matrix(marks.matrices[[i]]))
  
}
seriate.m <- seriate.m[complete.cases(seriate.m), ]


# 8. get not significant genes for each mark
ns.genes.list <- list()

for (i in 1:9) {
  
  ns.genes.list[[i]] <- my.function(mark = marks[i])
  
}


# 9. retrieve list of genes not significant for all marks
ns.genes <- ns.genes.list[[1]]

for (i in 2:9) {
  
  ns.genes <- intersect(ns.genes, ns.genes.list[[i]])
  
}

seriate.m <- seriate.m[!(rownames(seriate.m) %in% ns.genes), ]




# 10. seriation + heatmap
o <- seriate(seriate.m, method = "PCA_angle")
HT <-  Heatmap(seriate.m[, 1:12],
          row_order = get_order(o, 1),
          col = palette2, 
          name = marks[1], column_title = marks[1],
          column_title_gp = gpar(fontsize = 25),
          column_names_gp = gpar(fontsize = 25),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = T,
          heatmap_legend_param = list(title = "z-score"),
          gap = unit(5, "mm")) +
  
  Heatmap(seriate.m[, 13:24],
          row_order = get_order(o, 1),
          col = palette2, 
          name = marks[2], column_title = marks[2],
          column_title_gp = gpar(fontsize = 25),
          column_names_gp = gpar(fontsize = 25),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = T,
          heatmap_legend_param = list(title = "z-score"),
          gap = unit(5, "mm")) +
  
  Heatmap(seriate.m[, 25:36],
          row_order = get_order(o, 1),
          col = palette2, 
          name = marks[3], column_title = marks[3],
          column_title_gp = gpar(fontsize = 25),
          column_names_gp = gpar(fontsize = 25),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = T,
          heatmap_legend_param = list(title = "z-score"),
          gap = unit(5, "mm")) +
  
  Heatmap(seriate.m[, 37:48],
          row_order = get_order(o, 1),
          col = palette2, 
          name = marks[4], column_title = marks[4],
          column_title_gp = gpar(fontsize = 25),
          column_names_gp = gpar(fontsize = 25),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = T,
          heatmap_legend_param = list(title = "z-score"),
          gap = unit(5, "mm")) +
  
  Heatmap(seriate.m[, 49:60],
          row_order = get_order(o, 1),
          col = palette2, 
          name = marks[5], column_title = marks[5],
          column_title_gp = gpar(fontsize = 25),
          column_names_gp = gpar(fontsize = 25),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = T,
          heatmap_legend_param = list(title = "z-score"),
          gap = unit(5, "mm")) +
  
  Heatmap(seriate.m[, 61:72],
          row_order = get_order(o, 1),
          col = palette2, 
          name = marks[6], column_title = marks[6],
          column_title_gp = gpar(fontsize = 25),
          column_names_gp = gpar(fontsize = 25),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = T,
          heatmap_legend_param = list(title = "z-score"),
          gap = unit(5, "mm")) +
  
  Heatmap(seriate.m[, 73:84],
          row_order = get_order(o, 1),
          col = palette2, 
          name = marks[7], column_title = marks[7],
          column_title_gp = gpar(fontsize = 25),
          column_names_gp = gpar(fontsize = 25),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = T,
          heatmap_legend_param = list(title = "z-score"),
          gap = unit(5, "mm")) +
  
  Heatmap(seriate.m[, 85:96],
          row_order = get_order(o, 1),
          col = palette2, 
          name = marks[8], column_title = marks[8],
          column_title_gp = gpar(fontsize = 25),
          column_names_gp = gpar(fontsize = 25),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = T,
          heatmap_legend_param = list(title = "z-score"),
          gap = unit(5, "mm")) +
  
  Heatmap(seriate.m[, 97:108],
          row_order = get_order(o, 1),
          col = palette2, 
          name = marks[9], column_title = marks[9],
          column_title_gp = gpar(fontsize = 25),
          column_names_gp = gpar(fontsize = 25),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = T,
          heatmap_legend_param = list(title = "z-score"),
          gap = unit(5, "mm"))


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1t.pdf",
    width = 25, height = 3)
print(HT)
dev.off()


write.table(setdiff(silent.genes, ns.genes), 
            file = "~/public_html/Borsari_et_al_transdifferentiation_chromatin/genes/silent.genes.expression/silent.genes.variable.for.one.mark.txt",
            row.names = F, col.names = F, quote = F)

write.table(ns.genes, 
            file = "~/public_html/Borsari_et_al_transdifferentiation_chromatin/genes/silent.genes.expression/silent.genes.flat.for.all.marks.txt",
            row.names = F, col.names = F, quote = F)
