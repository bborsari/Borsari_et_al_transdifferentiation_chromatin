.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")



#************
# LIBRARIES *
#************


library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(circlize)


#************
# FUNCTIONS *
#************

my.function <- function(mark) {
        
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
        
        # 5. retrieve genes that are significantly variable for the mark
        mark.sig.genes <- read.table(paste0(mark, "/QN.merged/", mark, ".QN.merged.maSigPro.out.tsv"), h=T, sep="\t")
        mark.sig.genes$gene_id <- rownames(mark.sig.genes)
        mark.sig.genes <- mark.sig.genes %>% separate(gene_id, c("gene", "id"), "\\.")
        rownames(mark.sig.genes) <- mark.sig.genes$gene
        mark.sig.genes$gene <- NULL
        mark.sig.genes$id <- NULL
        mark.sig.genes <- rownames(mark.sig.genes)[rownames(mark.sig.genes) %in% rownames(expression.matrix)]
        
        # 6. retrieve genes that are flat for the mark
        mark.ns.genes <- setdiff(rownames(expression.matrix), mark.sig.genes)
        
        # 7. prepare separate expression and mark matrices for genes that are 
                # flat for expression only (sig)
                # flat for both expression and mark (ns)
        mark.matrix.scaled.sig <- mark.matrix.scaled[rownames(mark.matrix.scaled) %in% mark.sig.genes, ]
        expression.matrix.scaled.sig <- expression.matrix.scaled[rownames(expression.matrix.scaled) %in% mark.sig.genes, ]
        mark.matrix.scaled.ns <- mark.matrix.scaled[rownames(mark.matrix.scaled) %in% mark.ns.genes, ]
        expression.matrix.scaled.ns <- expression.matrix.scaled[rownames(expression.matrix.scaled) %in% mark.ns.genes, ]
        
        stopifnot(identical(rownames(mark.matrix.scaled.sig),
                            rownames(expression.matrix.scaled.sig)))
        
        stopifnot(identical(rownames(mark.matrix.scaled.ns),
                            rownames(expression.matrix.scaled.ns)))
        
        # 8. retrieve row order for sig. genes with clustering
        ht.sig <- Heatmap(mark.matrix.scaled.sig,
                          col = palette2, 
                          name = mark, column_title = mark,
                          column_title_gp = gpar(fontsize = 60),
                          column_names_gp = gpar(fontsize = 60),
                          show_row_names = FALSE, width = unit(60, "mm"),
                          cluster_rows = T, cluster_columns = F,
                          clustering_distance_rows = "euclidean",
                          clustering_method_rows = "complete",
                          show_column_names = F,
                          show_heatmap_legend = F,
                          heatmap_legend_param = list(title = "z-score"),
                          gap = unit(5, "mm"))
        
        new.order.sig <- rownames(mark.matrix.scaled.sig)[row_order(ht.sig)[[1]]]
        mark.matrix.scaled.sig <- mark.matrix.scaled.sig[new.order.sig, ]
        expression.matrix.scaled.sig <- expression.matrix.scaled.sig[new.order.sig, ]
        
        
        # 9. retrieve row order for ns genes with clustering
        # ht.ns <- Heatmap(mark.matrix.scaled.ns,
        #                  col = palette2, 
        #                  name = mark, column_title = mark,
        #                  column_title_gp = gpar(fontsize = 60),
        #                  column_names_gp = gpar(fontsize = 60),
        #                  show_row_names = FALSE, width = unit(60, "mm"),
        #                  cluster_rows = T, cluster_columns = F,
        #                  clustering_distance_rows = "euclidean",
        #                  clustering_method_rows = "complete",
        #                  show_column_names = F,
        #                  show_heatmap_legend = F,
        #                  heatmap_legend_param = list(title = "z-score"),
        #                  gap = unit(5, "mm"))
        
        # new.order.ns <- rownames(mark.matrix.scaled.ns)[row_order(ht.ns)[[1]]]
        new.order.ns <- rownames(mark.matrix.scaled.ns)
        mark.matrix.scaled.ns <- mark.matrix.scaled.ns[new.order.ns, ]
        expression.matrix.scaled.ns <- expression.matrix.scaled.ns[new.order.ns, ]
        
        
        # 10. plot the heatmap
        mark_ht_list = Heatmap(rbind(expression.matrix.scaled.sig, 
                                     expression.matrix.scaled.ns),
                               col = palette2,
                               name = "log2(TPM+1)", column_title = "expression",
                               column_title_gp = gpar(fontsize = 22),
                               column_names_gp = gpar(fontsize = 22),
                               show_row_names = FALSE, width = unit(40, "mm"),
                               cluster_rows = F, cluster_columns = F,
                               show_column_names = F,
                               show_heatmap_legend = F,
                               heatmap_legend_param = list(title = "z-score"),
                               split = c(rep("changes in mark", nrow(expression.matrix.scaled.sig)), 
                                         rep("no changes", nrow(expression.matrix.scaled.ns))),
                               gap = unit(5, "mm")) + #,
                               # border = T ) +
                
                Heatmap(rbind(mark.matrix.scaled.sig, 
                              mark.matrix.scaled.ns),
                        col = palette2, 
                        name = mark, column_title = mark,
                        column_title_gp = gpar(fontsize = 22),
                        column_names_gp = gpar(fontsize = 22),
                        show_row_names = FALSE, width = unit(40, "mm"),
                        cluster_rows = F, cluster_columns = F,
                        show_column_names = F,
                        show_heatmap_legend = F,
                        heatmap_legend_param = list(title = "z-score"),
                        split = c(rep("changes in mark", nrow(expression.matrix.scaled.sig)), 
                                  rep("no changes", nrow(expression.matrix.scaled.ns))),
                        gap = unit(5, "mm")) #,
                        # border = T)

        
        return(list(mark_ht_list, new.order.sig, new.order.ns))
        
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

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1n.H3K27ac.pdf", height=5, width=6)
H3K27ac <- my.function(mark = "H3K27ac")

# write.table(as.data.frame(H3K27ac[[2]]), 
#             "~/public_html/Borsari_et_al_transdifferentiation_chromatin/genes/flat.genes.expression/H3K27ac.txt",
#             col.names = F, row.names = F, quote=F, sep="\t")
# write.table(as.data.frame(H3K27ac[[3]]), 
#             "~/public_html/Borsari_et_al_transdifferentiation_chromatin/genes/flat.genes.expression/H3K27ac.2.txt",
#             col.names = F, row.names = F, quote=F, sep="\t")


print(H3K27ac)
dev.off()


###--------
### H3K9ac
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1n.H3K9ac.pdf", height=5, width=6)
H3K9ac <- my.function(mark="H3K9ac")
print(H3K9ac)
dev.off()


###--------
### H4K20me1
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1n.H4K20me1.pdf", height=5, width=6)
H4K20me1 <- my.function(mark="H4K20me1")
print(H4K20me1)
dev.off()


###--------
### H3K4me3
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1n.H3K4me3.pdf", height=5, width=6)
H3K4me3 <- my.function(mark = "H3K4me3")
print(H3K4me3)
dev.off()


###--------
### H3K4me1
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1n.H3K4me1.pdf", height=5, width=6)
H3K4me1 <- my.function(mark="H3K4me1")
print(H3K4me1)
dev.off()


###--------
### H3K36me3
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1n.H3K36me3.pdf", height=5, width=6)
H3K36me3 <- my.function(mark="H3K36me3")
print(H3K36me3)
dev.off()


###--------
### H3K4me2
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1n.H3K4me2.pdf", height=5, width=6)
H3K4me2 <- my.function(mark="H3K4me2")
print(H3K4me2)
dev.off()


###--------
### H3K9me3
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1n.H3K9me3.pdf", height=5, width=6)
H3K9me3 <- my.function(mark="H3K9me3")
print(H3K9me3)
dev.off()


###--------
### H3K27me3
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1n.H3K27me3.pdf", height=5, width=6)
H3K27me3 <- my.function(mark="H3K27me3")
print(H3K27me3)
dev.off()
