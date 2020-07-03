


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
  
  # 1. read mark matrix
  mark.matrix <- read.table(paste0(mark, "/QN.merged/", mark, ".matrix.tsv"), h=T, sep="\t")
  
  # 2. reorder mark matrix according to metadata
  mark.matrix <- mark.matrix[rownames(metadata), ]
  stopifnot(identical(rownames(metadata), rownames(mark.matrix)))
  
  # 3. center and scale mark matrix
  mark.matrix.scaled <- as.data.frame(t(scale(t(mark.matrix))))
  
  # 4. read matrix of 6 groups
  mark.6.groups <- read.table(paste0(mark, "/QN.merged/", mark, ".6.groups.tsv"),
                                 h=T, sep="\t")
  mark.6.groups$group <- gsub("peak_not_TSS", "no_peak", mark.6.groups$group)
  
  
  # 5. genes w/o peaks of the mark
  mark.no.peak.genes <- rownames(mark.6.groups[mark.6.groups$group == "no_peak", ])
  mark.metadata.no.peak.genes <- metadata[rownames(metadata) %in% mark.no.peak.genes, ]
  stopifnot(identical(rownames(expression.matrix.scaled[rownames(expression.matrix.scaled) %in% rownames(mark.metadata.no.peak.genes), ]),
                      rownames(mark.matrix.scaled[rownames(mark.matrix.scaled) %in% rownames(mark.metadata.no.peak.genes), ])))
  
  
  
  # 6. genes w/ peak of the mark and stable profile
  mark.ns.genes <- rownames(mark.6.groups[mark.6.groups$group == "stable", ])
  mark.metadata.ns.genes <- metadata[rownames(metadata) %in% mark.ns.genes, ]
  stopifnot(identical(rownames(expression.matrix.scaled[rownames(expression.matrix.scaled) %in% rownames(mark.metadata.ns.genes), ]),
                      rownames(mark.matrix.scaled[rownames(mark.matrix.scaled) %in% rownames(mark.metadata.ns.genes), ])))
  
  
  # 6. genes w/ peak of the mark and variable profile
  mark.sig.genes <- rownames(mark.6.groups[mark.6.groups$group %in% c("positively_correlated",
                                                                      "negatively_correlated",
                                                                      "not_correlated"), ])
  mark.metadata.sig.genes <- metadata[rownames(metadata) %in% mark.sig.genes, ]
  stopifnot(identical(rownames(expression.matrix.scaled[rownames(expression.matrix.scaled) %in% rownames(mark.metadata.sig.genes), ]),
                      rownames(mark.matrix.scaled[rownames(mark.matrix.scaled) %in% rownames(mark.metadata.sig.genes), ])))
  
  mark.partition1 = rbind(mark.metadata.sig.genes,
                          mark.metadata.ns.genes,
                          mark.metadata.no.peak.genes)[, "final_class"]
  mark.partition2 = rbind(mark.metadata.sig.genes,
                          mark.metadata.ns.genes,
                          mark.metadata.no.peak.genes)[, "avg_exp"]
  
  colnames(expression.matrix.scaled) <- paste0(c(0, 3, 6, 9, 12, 18, 24, 36, 48, 72, 120, 168), "h")
  colnames(mark.matrix.scaled) <- paste0(c(0, 3, 6, 9, 12, 18, 24, 36, 48, 72, 120, 168), "h")
  
  
  
  mark.ht_list = Heatmap(mark.partition1, 
                         name = "profiles", 
                         col= palette[names(palette) %in% c("bending",
                                                            "down-regulated",
                                                            "peaking",
                                                            "up-regulated")],
                         show_heatmap_legend = F,
                         show_row_names = FALSE, width = unit(15, "mm"),
                         show_column_names = F,
                         split = c(rep("changes in \nexpression & mark", nrow(mark.metadata.sig.genes)), 
                                   rep("changes in \nexpression only", nrow(mark.metadata.ns.genes)),
                                   rep("no_peak", nrow(mark.metadata.no.peak.genes))),
                         gap = unit(4, "mm")) +

    Heatmap(rbind(expression.matrix.scaled[rownames(expression.matrix.scaled) %in% rownames(mark.metadata.sig.genes), ], 
                  expression.matrix.scaled[rownames(expression.matrix.scaled) %in% rownames(mark.metadata.ns.genes), ],
                  expression.matrix.scaled[rownames(expression.matrix.scaled) %in% rownames(mark.metadata.no.peak.genes), ]),
            col = palette2,
            name = "log2(TPM+1)", column_title = "expression",
            column_title_gp = gpar(fontsize = 60),
            column_names_gp = gpar(fontsize = 60),
            show_row_names = FALSE, width = unit(120, "mm"),
            cluster_rows = F, cluster_columns = F,
            show_column_names = F,
            show_heatmap_legend = F,
            heatmap_legend_param = list(title = "z-score"),
            split = c(rep("changes in \nexpression & mark", nrow(mark.metadata.sig.genes)), 
                      rep("changes in \nexpression only", nrow(mark.metadata.ns.genes)),
                      rep("no_peak", nrow(mark.metadata.no.peak.genes))),
            gap = unit(4, "mm"),
            border = T) +
    
    Heatmap(rbind(mark.matrix.scaled[rownames(mark.matrix.scaled) %in% rownames(mark.metadata.sig.genes), ], 
                  mark.matrix.scaled[rownames(mark.matrix.scaled) %in% rownames(mark.metadata.ns.genes), ],
                  mark.matrix.scaled[rownames(mark.matrix.scaled) %in% rownames(mark.metadata.no.peak.genes), ]),
            col = palette2, 
            name = mark, column_title = mark,
            column_title_gp = gpar(fontsize = 60),
            column_names_gp = gpar(fontsize = 60),
            show_row_names = FALSE, width = unit(120, "mm"),
            cluster_rows = F, cluster_columns = F,
            show_column_names = F,
            show_heatmap_legend = F,
            heatmap_legend_param = list(title = "z-score"),
            split = c(rep("changes in \nexpression & mark", nrow(mark.metadata.sig.genes)), 
                      rep("changes in \nexpression only", nrow(mark.metadata.ns.genes)),
                      rep("no_peak", nrow(mark.metadata.no.peak.genes))),
            gap = unit(4, "mm"),
            border = T)
  
  return(mark.ht_list)
  
}


#********
# BEGIN *
#********

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/")

palette <- c("down-regulated" = "#2d7f89", 
             "bending" = "#7acbd5",
             "up-regulated" = "#89372d",
             "peaking" = "#d5847a"
)

palette2 = rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061'))


# 1. import expression matrix
expression.matrix <- read.table("expression/QN.merged/expression.matrix.tsv", h=T, sep="\t")

# 2. import metadata and check row order is the same as in expression matrix
metadata <- read.table("expression/QN.merged/metadata.class2.tsv", h=T, sep="\t")
stopifnot(identical(rownames(expression.matrix), rownames(metadata)))

# 3. work with final_class column
metadata$final_class <- gsub("regulation", "-regulated", metadata$final_class)

# 4. add to metadata info about average expression and time-point
metadata$avg_exp <- apply(expression.matrix, 1, mean)
metadata$tp1 <- apply(expression.matrix, 1, which.max)
metadata$tp2 <- apply(expression.matrix, 1, which.min)
metadata$tp <- ifelse(metadata$final_class %in% c("bending", "down-regulated"), 
                      metadata$tp2, 
                      metadata$tp1)
metadata <- metadata[order(metadata$final_class, metadata$tp, metadata$avg_exp), ]

# 5. reorder expression matrix according to metadata row order
expression.matrix <- expression.matrix[rownames(metadata), ]
stopifnot(identical(rownames(expression.matrix), rownames(metadata)))

# 6. center and scale expression matrix
expression.matrix.scaled <- as.data.frame(t(scale(t(expression.matrix))))


###----------
### H3K27ac
###----------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1l2.H3K27ac.pdf", height=16, width=12)
H3K27ac <- my.function(mark = "H3K27ac")
print(H3K27ac)
dev.off()


###--------
### H3K9ac
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1l2.H3K9ac.pdf", height=16, width=12)
H3K9ac <- my.function(mark="H3K9ac")
print(H3K9ac)
dev.off()


###--------
### H4K20me1
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1l2.H4K20me1.pdf", height=16, width=12)
H4K20me1 <- my.function(mark="H4K20me1")
print(H4K20me1)
dev.off()


###--------
### H3K4me3
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1l2.H3K4me3.pdf", height=16, width=12)
H3K4me3 <- my.function(mark = "H3K4me3")
print(H3K4me3)
dev.off()


###--------
### H3K4me1
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1l2.H3K4me1.pdf", height=16, width=12)
H3K4me1 <- my.function(mark="H3K4me1")
print(H3K4me1)
dev.off()


###--------
### H3K36me3
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1l2.H3K36me3.pdf", height=16, width=12)
H3K36me3 <- my.function(mark="H3K36me3")
print(H3K36me3)
dev.off()


###--------
### H3K4me2
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1l2.H3K4me2.pdf", height=16, width=12)
H3K4me2 <- my.function(mark="H3K4me2")
print(H3K4me2)
dev.off()


###--------
### H3K9me3
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1l2.H3K9me3.pdf", height=16, width=12)
H3K9me3 <- my.function(mark="H3K9me3")
print(H3K9me3)
dev.off()


###--------
### H3K27me3
###--------

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1l2.H3K27me3.pdf", height=16, width=12)
H3K27me3 <- my.function(mark="H3K27me3")
print(H3K27me3)
dev.off()
