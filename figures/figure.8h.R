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
  
  # 1. read mark matrix of all genes after QN merged normalization
  mark.matrix <- read.table(paste0(mark, "/QN.merged/", mark, ".matrix.after.QN.merged.tsv"),
                            h=T, sep="\t")
  
  # 2. keep only flat genes (aka genes in expression matrix)
  mark.matrix <- mark.matrix[rownames(mark.matrix) %in% rownames(expression.matrix.scaled), ]
  
  # 3. check order of rows is the same as with expression
  mark.matrix <- mark.matrix[rownames(expression.matrix), ]
  stopifnot(identical(rownames(expression.matrix.scaled), rownames(mark.matrix)))
  
  # 4. center and scale mark matrix of flat genes 
  mark.matrix.scaled <- as.data.frame(t(scale(t(mark.matrix))))
  
  # 5. genes w/ peaks of the mark
  mark.genes.peaks <- read.table(paste0(mark, "/QN.merged/all.genes.intersecting.peaks.tsv"),
                                 h=F, sep="\t", stringsAsFactors = F)
  mark.genes.peaks <- mark.genes.peaks$V1
  
  # 6. genes w/ variable profile for the mark
  mark.sig.genes <- read.table(paste0(mark, "/QN.merged/", mark, ".QN.merged.maSigPro.out.tsv"), 
                               h=T, sep="\t")
  mark.sig.genes$gene_id <- rownames(mark.sig.genes)
  mark.sig.genes <- mark.sig.genes %>% separate(gene_id, c("gene", "id"), "\\.")
  rownames(mark.sig.genes) <- mark.sig.genes$gene
  mark.sig.genes$gene <- NULL
  mark.sig.genes$id <- NULL
  mark.sig.genes <- rownames(mark.sig.genes)[rownames(mark.sig.genes) %in% rownames(expression.matrix.scaled)]
  
  # 7. genes w/ peak of the mark and variable profile
  mark.sig.genes <- intersect(mark.sig.genes, mark.genes.peaks)
  stopifnot(identical(rownames(expression.matrix.scaled[rownames(expression.matrix.scaled) %in% mark.sig.genes, ]),
                      rownames(mark.matrix.scaled[rownames(mark.matrix.scaled) %in% mark.sig.genes, ])))
  
  print(mark)
  print(length(mark.sig.genes))
  
  
  ## 7.1. apply seriation on mark matrix only
  o.b <- seriate(as.matrix(mark.matrix.scaled[rownames(mark.matrix.scaled) %in% mark.sig.genes, ]),
                 method = "PCA_angle")
  
  ## 7.2. reorder corresponding expression and mark matrix according to seriation order
  exp.b <- expression.matrix.scaled[rownames(expression.matrix.scaled) %in% mark.sig.genes, ]
  exp.b <- exp.b[get_order(o.b, 1), ]
  
  mark.b <- mark.matrix.scaled[rownames(mark.matrix.scaled) %in% mark.sig.genes, ]
  mark.b <- mark.b[get_order(o.b, 1), ]
  stopifnot(identical(rownames(exp.b), rownames(mark.b)))
  
  
  # 8. heatmap
  colnames(exp.b) <- paste0(c(0, 3, 6, 9, 12, 18, 24, 36, 48, 72, 120, 168), "h")
  colnames(mark.b) <- paste0(c(0, 3, 6, 9, 12, 18, 24, 36, 48, 72, 120, 168), "h")
  
  
  mark.ht_list = 
    
    Heatmap(exp.b,
            col = palette2,
            name = "log2(TPM+1)", column_title = "expression",
            column_title_gp = gpar(fontsize = 60),
            column_names_gp = gpar(fontsize = 60),
            show_row_names = FALSE, width = unit(120, "mm"),
            cluster_rows = F, cluster_columns = F,
            show_column_names = F,
            show_heatmap_legend = F,
            heatmap_legend_param = list(title = "z-score")) +
    
    Heatmap(mark.b,
            col = palette2, 
            name = mark, column_title = mark,
            column_title_gp = gpar(fontsize = 60),
            column_names_gp = gpar(fontsize = 60),
            show_row_names = FALSE, width = unit(120, "mm"),
            cluster_rows = F, cluster_columns = F,
            show_column_names = F,
            show_heatmap_legend = F,
            heatmap_legend_param = list(title = "z-score"))
  
  return(mark.ht_list)
  
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


###----------
### H3K27ac
###----------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8h.H3K27ac.pdf", height=6, width=12)
H3K27ac <- my.function(mark = "H3K27ac")
print(H3K27ac)
dev.off()


###--------
### H3K9ac
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8h.H3K9ac.pdf", height=6, width=12)
H3K9ac <- my.function(mark="H3K9ac")
print(H3K9ac)
dev.off()


###--------
### H4K20me1
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8h.H4K20me1.pdf", height=6, width=12)
H4K20me1 <- my.function(mark="H4K20me1")
print(H4K20me1)
dev.off()


###--------
### H3K4me3
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8h.H3K4me3.pdf", height=6, width=12)
H3K4me3 <- my.function(mark = "H3K4me3")
print(H3K4me3)
dev.off()


###--------
### H3K4me1
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8h.H3K4me1.pdf", height=6, width=12)
H3K4me1 <- my.function(mark="H3K4me1")
print(H3K4me1)
dev.off()


###--------
### H3K36me3
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8h.H3K36me3.pdf", height=6, width=12)
H3K36me3 <- my.function(mark="H3K36me3")
print(H3K36me3)
dev.off()


###--------
### H3K4me2
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8h.H3K4me2.pdf", height=6, width=12)
H3K4me2 <- my.function(mark="H3K4me2")
print(H3K4me2)
dev.off()


###--------
### H3K9me3
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8h.H3K9me3.pdf", height=6, width=12)
H3K9me3 <- my.function(mark="H3K9me3")
print(H3K9me3)
dev.off()


###--------
### H3K27me3
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8h.H3K27me3.pdf", height=6, width=12)
H3K27me3 <- my.function(mark="H3K27me3")
print(H3K27me3)
dev.off()
