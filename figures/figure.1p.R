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




my.function2 <- function(mark) {
  
  # 1. retrieve genes that are significantly variable for the mark
  # while flat for expression
  mark.sig.genes <- read.table(paste0(mark, "/QN.merged/", mark, ".QN.merged.maSigPro.out.tsv"), 
                               h=T, sep="\t")
  mark.sig.genes$gene_id <- rownames(mark.sig.genes)
  mark.sig.genes <- mark.sig.genes %>% separate(gene_id, c("gene", "id"), "\\.")
  rownames(mark.sig.genes) <- mark.sig.genes$gene
  mark.sig.genes$gene <- NULL
  mark.sig.genes$id <- NULL
  mark.sig.genes <- rownames(mark.sig.genes)[rownames(mark.sig.genes) %in% rownames(expression.matrix)]
  
  # 2. for each mark keep only the selected variable genes 
  tmp.list <- list()
  for ( i in 1:10 ){
    
    tmp <- marks.matrices[[i]]
    tmp <- tmp[rownames(tmp) %in% mark.sig.genes, ]
    tmp.list[[i]] <- tmp
    
  }
  
  # 3. build input matrix for seriation only with the selected genes
  seriate.m <- as.matrix(tmp.list[[1]])
  
  for ( i in 2:10 ) {
    
    seriate.m <- cbind(seriate.m, as.matrix(tmp.list[[i]]))
    
  }
  
  seriate.m <- seriate.m[complete.cases(seriate.m), ]
  
  return(seriate.m)
  

}


my.function3 <- function(f3.method, seriate.m) {
  
  # 1. seriation (without expression)
  o <- seriate(seriate.m[, 13:ncol(seriate.m)], method = f3.method)
  
  # 2. heatmap
  HT <-  Heatmap(seriate.m[, 1:12],
                 row_order = get_order(o, 1),
                 col = palette2, 
                 name = "expression", column_title = "expression",
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
            name = marks[1], column_title = marks[1],
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
            name = marks[2], column_title = marks[2],
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
            name = marks[3], column_title = marks[3],
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
            name = marks[4], column_title = marks[4],
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
            name = marks[5], column_title = marks[5],
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
            name = marks[6], column_title = marks[6],
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
            name = marks[7], column_title = marks[7],
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
            name = marks[8], column_title = marks[8],
            column_title_gp = gpar(fontsize = 25),
            column_names_gp = gpar(fontsize = 25),
            show_row_names = FALSE, width = unit(50, "mm"),
            cluster_rows = F, cluster_columns = F,
            show_column_names = F,
            show_heatmap_legend = T,
            heatmap_legend_param = list(title = "z-score"),
            gap = unit(5, "mm")) +
    
    Heatmap(seriate.m[, 109:120],
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
  
    
  return(HT)
  
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
marks <- c("H3K4me1", "H3K4me2", "H3K4me3", 
           "H3K9ac", "H3K27ac", "H4K20me1", 
           "H3K36me3", "H3K27me3", "H3K9me3")
marks.matrices <- list()

marks.matrices[[1]] <- expression.matrix.scaled

for ( i in 1:9 ) {
  
  marks.matrices[[i+1]] <- my.function1(mark = marks[i])
  
}


# 7. get input matrix for seriation 
# selecting sig. genes for each mark
list_seriation_methods("matrix")


#----------
# H3K27ac
#----------

H3K27ac.seriate <- my.function2(mark="H3K27ac")

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1p.H3K27ac.PCA.pdf",
    width=25, height=6)
print(my.function3(f3.method = "PCA", seriate.m = H3K27ac.seriate))
dev.off()
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1p.H3K27ac.PCA_angle.pdf",
    width=25, height=6)
print(my.function3(f3.method = "PCA_angle", seriate.m = H3K27ac.seriate))
dev.off()



#----------
# H3K9ac
#----------

H3K9ac.seriate <- my.function2(mark="H3K9ac")

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1p.H3K9ac.PCA.pdf",
    width=25, height=6)
print(my.function3(f3.method = "PCA", seriate.m = H3K9ac.seriate))
dev.off()
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1p.H3K9ac.PCA_angle.pdf",
    width=25, height=6)
print(my.function3(f3.method = "PCA_angle", seriate.m = H3K9ac.seriate))
dev.off()




#----------
# H4K20me1
#----------

H4K20me1.seriate <- my.function2(mark="H4K20me1")

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1p.H4K20me1.PCA.pdf",
    width=25, height=6)
print(my.function3(f3.method = "PCA", seriate.m = H4K20me1.seriate))
dev.off()

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1p.H4K20me1.PCA_angle.pdf",
    width=25, height=6)
print(my.function3(f3.method = "PCA_angle", seriate.m = H4K20me1.seriate))
dev.off()



#----------
# H3K4me1
#----------

H3K4me1.seriate <- my.function2(mark="H3K4me1")

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1p.H3K4me1.PCA.pdf",
    width=25, height=6)
print(my.function3(f3.method = "PCA", seriate.m = H3K4me1.seriate))
dev.off()

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1p.H3K4me1.PCA_angle.pdf",
    width=25, height=6)
print(my.function3(f3.method = "PCA_angle", seriate.m = H3K4me1.seriate))
dev.off()



#----------
# H3K4me3
#----------

H3K4me3.seriate <- my.function2(mark="H3K4me3")

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1p.H3K4me3.PCA.pdf",
    width=25, height=6)
print(my.function3(f3.method = "PCA", seriate.m = H3K4me3.seriate))
dev.off()

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1p.H3K4me3.PCA_angle.pdf",
    width=25, height=6)
print(my.function3(f3.method = "PCA_angle", seriate.m = H3K4me3.seriate))
dev.off()



#----------
# H3K4me2
#----------

H3K4me2.seriate <- my.function2(mark="H3K4me2")

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1p.H3K4me2.PCA.pdf",
    width=25, height=6)
print(my.function3(f3.method = "PCA", seriate.m = H3K4me2.seriate))
dev.off()

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1p.H3K4me2.PCA_angle.pdf",
    width=25, height=6)
print(my.function3(f3.method = "PCA_angle", seriate.m = H3K4me2.seriate))
dev.off()



#----------
# H3K9me3
#----------

H3K9me3.seriate <- my.function2(mark="H3K9me3")

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1p.H3K9me3.PCA.pdf",
    width=25, height=6)
print(my.function3(f3.method = "PCA", seriate.m = H3K9me3.seriate))
dev.off()

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1p.H3K9me3.PCA_angle.pdf",
    width=25, height=6)
print(my.function3(f3.method = "PCA_angle", seriate.m = H3K9me3.seriate))
dev.off()



#----------
# H3K36me3
#----------

H3K36me3.seriate <- my.function2(mark="H3K36me3")

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1p.H3K36me3.PCA.pdf",
    width=25, height=6)
print(my.function3(f3.method = "PCA", seriate.m = H3K36me3.seriate))
dev.off()

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1p.H3K36me3.PCA_angle.pdf",
    width=25, height=6)
print(my.function3(f3.method = "PCA_angle", seriate.m = H3K36me3.seriate))
dev.off()



#----------
# H3K27me3
#----------

H3K27me3.seriate <- my.function2(mark="H3K27me3")

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1p.H3K27me3.PCA.pdf",
    width=25, height=6)
print(my.function3(f3.method = "PCA", seriate.m = H3K27me3.seriate))
dev.off()

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1p.H3K27me3.PCA_angle.pdf",
    width=25, height=6)
print(my.function3(f3.method = "PCA_angle", seriate.m = H3K27me3.seriate))
dev.off()