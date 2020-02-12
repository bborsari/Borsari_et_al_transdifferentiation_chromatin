


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
  
  mark.matrix <- read.table(paste0(mark, "/QN.merged/", mark, ".matrix.tsv"), h=T, sep="\t")
  mark.matrix <- mark.matrix[rownames(metadata), ]
  stopifnot(identical(rownames(metadata), rownames(mark.matrix)))
  mark.matrix.scaled <- as.data.frame(t(scale(t(mark.matrix))))
  mark.sig.genes <- read.table(paste0(mark, "/QN.merged/", mark, ".QN.merged.maSigPro.out.tsv"), h=T, sep="\t")
  mark.sig.genes$gene_id <- rownames(mark.sig.genes)
  mark.sig.genes <- mark.sig.genes %>% separate(gene_id, c("gene", "id"), "\\.")
  rownames(mark.sig.genes) <- mark.sig.genes$gene
  mark.sig.genes$gene <- NULL
  mark.sig.genes$id <- NULL
  mark.sig.genes <- mark.sig.genes[rownames(mark.sig.genes) %in% rownames(expression.matrix), , drop=F]
  
  mark.metadata.sig.genes <- metadata[rownames(metadata) %in% rownames(mark.sig.genes), ]
  mark.ns.genes <- setdiff(rownames(expression.matrix), rownames(mark.sig.genes))
  mark.metadata.ns.genes <- metadata[rownames(metadata) %in% mark.ns.genes, ]
  
  stopifnot(identical(rownames(expression.matrix.scaled[rownames(expression.matrix.scaled) %in% rownames(mark.metadata.sig.genes), ]),
                      rownames(mark.matrix.scaled[rownames(mark.matrix.scaled) %in% rownames(mark.metadata.sig.genes), ])))
  stopifnot(identical(rownames(expression.matrix.scaled[rownames(expression.matrix.scaled) %in% rownames(mark.metadata.ns.genes), ]),
                      rownames(mark.matrix.scaled[rownames(mark.matrix.scaled) %in% rownames(mark.metadata.ns.genes), ])))
  
  mark.partition1 = rbind(mark.metadata.sig.genes,
                          mark.metadata.ns.genes)[, "class4"]
  mark.partition2 = rbind(mark.metadata.sig.genes,
                          mark.metadata.ns.genes)[, "avg_exp"]
  
  colnames(expression.matrix.scaled) <- paste0(c(0, 3, 6, 9, 12, 18, 24, 36, 48, 72, 120, 168), "h")
  colnames(mark.matrix.scaled) <- paste0(c(0, 3, 6, 9, 12, 18, 24, 36, 48, 72, 120, 168), "h")
  
  
  
  mark.ht_list = Heatmap(mark.partition1, 
                         name = "profiles", 
                         col= palette[names(palette) %in% c("bending",
                                                            "down-regulated",
                                                            "peaking",
                                                            "up-regulated")],
                         show_heatmap_legend = T,
                         show_row_names = FALSE, width = unit(15, "mm"),
                         show_column_names = F,
                         split = c(rep("changes in \nexpression & mark", nrow(mark.metadata.sig.genes)), 
                                   rep("changes in \nexpression only", nrow(mark.metadata.ns.genes))),
                         gap = unit(5, "mm")) +
    # Heatmap(mark.partition2,
    #         name = "mean exp. mark",
    #         cluster_rows = F,
    #         show_heatmap_legend = F,
    #         show_row_names = FALSE, width = unit(6, "mm"),
    #         show_column_names = F,
    #         split = c(rep("changes in \nexpression & mark", nrow(mark.metadata.sig.genes)),
    #                   rep("changes in \nexpression only", nrow(mark.metadata.ns.genes))),
    #         gap = unit(5, "mm"),) +
    
    Heatmap(rbind(expression.matrix.scaled[rownames(expression.matrix.scaled) %in% rownames(mark.metadata.sig.genes), ], 
                  expression.matrix.scaled[rownames(expression.matrix.scaled) %in% rownames(mark.metadata.ns.genes), ]),
            col = palette2,
            name = "log2(TPM+1)", column_title = "expression",
            column_title_gp = gpar(fontsize = 60),
            column_names_gp = gpar(fontsize = 60),
            show_row_names = FALSE, width = unit(50, "mm"),
            cluster_rows = F, cluster_columns = F,
            show_column_names = F,
            show_heatmap_legend = T,
            heatmap_legend_param = list(title = "z-score"),
            # bottom_annotation = HeatmapAnnotation(
            #   text = anno_text(colnames(expression.matrix.scaled), 
            #                    rot = 60, location = unit(1, "npc"), 
            #                    just = "right", gp = gpar(fontsize = 30)),
            #   annotation_height = max_text_width(colnames(expression.matrix.scaled))),
            split = c(rep("changes in \nexpression & mark", nrow(mark.metadata.sig.genes)), 
                      rep("changes in \nexpression only", nrow(mark.metadata.ns.genes))),
            gap = unit(5, "mm"),
            border = T) +
    Heatmap(rbind(mark.matrix.scaled[rownames(mark.matrix.scaled) %in% rownames(mark.metadata.sig.genes), ], 
                  mark.matrix.scaled[rownames(mark.matrix.scaled) %in% rownames(mark.metadata.ns.genes), ]),
            col = palette2, 
            name = mark, column_title = mark,
            column_title_gp = gpar(fontsize = 60),
            column_names_gp = gpar(fontsize = 60),
            show_row_names = FALSE, width = unit(50, "mm"),
            cluster_rows = F, cluster_columns = F,
            show_column_names = F,
            show_heatmap_legend = F,
            heatmap_legend_param = list(title = "z-score"),
            # bottom_annotation = HeatmapAnnotation(
            #   text = anno_text(colnames(expression.matrix.scaled), 
            #                    rot = 60, location = unit(1, "npc"), 
            #                    just = "right", gp = gpar(fontsize = 30)),
            #   annotation_height = max_text_width(colnames(expression.matrix.scaled))),
            split = c(rep("changes in \nexpression & mark", nrow(mark.metadata.sig.genes)), 
                      rep("changes in \nexpression only", nrow(mark.metadata.ns.genes))),
            gap = unit(5, "mm"),
            border = T)
  
  return(mark.ht_list)
  
}


#********
# BEGIN *
#********

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/")

palette <- c("down-regulated" = "#810f7c", 
             "bending" = "#a4a4a4",
             "up-regulated" = "#f7a673",
             "peaking" = "#c7e9b4"
)

palette2 = rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061'))


## 1. import expression and metadata matrices
expression.matrix <- read.table("expression/QN.merged/expression.matrix.tsv", h=T, sep="\t")
metadata <- read.table("expression/QN.merged/metadata.class2.tsv", h=T, sep="\t")
stopifnot(identical(rownames(expression.matrix), rownames(metadata)))
metadata$class <- gsub("regulation", "-regulated", metadata$class)
metadata$class2 <- gsub("regulation", "-regulated", metadata$class2)
metadata$class3 <- ifelse(metadata$class %in% c("down-regulated", "up-regulated", "peaking", "bending"),
                          as.character(metadata$class), paste0(metadata$class2, "#2"))
metadata$class4 <- ifelse(metadata$class %in% c("down-regulated", "up-regulated", "peaking", "bending"),
                          as.character(metadata$class), as.character(metadata$class2))

metadata$avg_exp <- apply(expression.matrix, 1, mean)
metadata$tp1 <- apply(expression.matrix, 1, which.max)
metadata$tp2 <- apply(expression.matrix, 1, which.min)
metadata$tp <- ifelse(metadata$class4 %in% c("bending", "down-regulated"), metadata$tp2, metadata$tp1)
metadata <- metadata[order(metadata$class4, metadata$tp, metadata$avg_exp), ]

expression.matrix <- expression.matrix[rownames(metadata), ]
expression.matrix.scaled <- as.data.frame(t(scale(t(expression.matrix))))


###----------
### H3K27ac
###----------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1l.H3K27ac.legend.pdf", height=15, width=12)
H3K27ac <- my.function(mark = "H3K27ac")
print(H3K27ac)
dev.off()


# ###--------
# ### H3K9ac
# ###--------
# 
# pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1l.H3K9ac.pdf", height=15, width=12)
# H3K9ac <- my.function(mark="H3K9ac")
# print(H3K9ac)
# dev.off()
# 
# 
# ###--------
# ### H4K20me1
# ###--------
# 
# pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1l.H4K20me1.pdf", height=15, width=12)
# H4K20me1 <- my.function(mark="H4K20me1")
# print(H4K20me1)
# dev.off()
# 
# 
# ###--------
# ### H3K4me3
# ###--------
# 
# pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1l.H3K4me3.pdf", height=15, width=12)
# H3K4me3 <- my.function(mark = "H3K4me3")
# print(H3K4me3)
# dev.off()
# 
# 
# ###--------
# ### H3K4me1
# ###--------
# 
# pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1l.H3K4me1.pdf", height=15, width=12)
# H3K4me1 <- my.function(mark="H3K4me1")
# print(H3K4me1)
# dev.off()
# 
# 
# ###--------
# ### H3K36me3
# ###--------
# 
# pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1l.H3K36me3.pdf", height=15, width=12)
# H3K36me3 <- my.function(mark="H3K36me3")
# print(H3K36me3)
# dev.off()
# 
# 
# ###--------
# ### H3K4me2
# ###--------
# 
# pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1l.H3K4me2.pdf", height=15, width=12)
# H3K4me2 <- my.function(mark="H3K4me2")
# print(H3K4me2)
# dev.off()
# 
# 
# ###--------
# ### H3K9me3
# ###--------
# 
# pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1l.H3K9me3.pdf", height=15, width=12)
# H3K9me3 <- my.function(mark="H3K9me3")
# print(H3K9me3)
# dev.off()
# 
# 
# ###--------
# ### H3K27me3
# ###--------
# 
# pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1l.H3K27me3.pdf", height=15, width=12)
# H3K27me3 <- my.function(mark="H3K27me3")
# print(H3K27me3)
# dev.off()
