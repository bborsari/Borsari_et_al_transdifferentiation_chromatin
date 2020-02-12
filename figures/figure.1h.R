.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


# plots of expression profiles of 4 example genes: peaking, upregulated, downregulated, bending

#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(ggpubr)
library(reshape2)
library(facetscales)

palette <- c("down-\nregulated\n\n(IGLL1)" = "#810f7c", 
             "bending\n\n(RGS1)" = "#737373",
             "up-\nregulated\n\n(C1QA)" = "#f16913",
             "peaking\n\n(RFX8)" = "#c7e9b4",
             "flat\n\n(ARF5)" = "#1c9099"
)

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/expression/")


#********
# BEGIN *
#********

expression.matrix <- read.table("QN.merged/selected.genes.rep.2.3.after.QN.merged.tsv", h=T, sep="\t")
metadata <- read.table("QN.merged/metadata.tsv", h=T, sep="\t")
hours <- c(0, 3, 6, 9, 12, 18, 24, 36, 48, 72, 120, 168)

# 1. compute FC
expression.matrix$FC <- apply(expression.matrix, 1, function(x){max(x) - min(x)})
expression.matrix <- expression.matrix[order(expression.matrix$FC, decreasing = T), ]

# 2. select upregulated gene - ENSG00000173372
upregulated.gene <- head(expression.matrix[rownames(expression.matrix) %in% rownames(metadata[metadata$class == "upregulation", ]), ], )
upregulated.gene <- tail(upregulated.gene, 1)
upregulated.gene <- melt(upregulated.gene[, 1:12])
upregulated.gene$tp <- hours
upregulated.gene$class <- "up-\nregulated\n\n(C1QA)"

# 3. select downregulated gene - ENSG00000128322
downregulated.gene <- head(expression.matrix[rownames(expression.matrix) %in% rownames(metadata[metadata$class == "downregulation", ]), ], 1)
downregulated.gene <- melt(downregulated.gene[, 1:12])
downregulated.gene$tp <- hours
downregulated.gene$class <- "down-\nregulated\n\n(IGLL1)"

# 4. select peaking gene
# peaking.gene <- expression.matrix["ENSG00000079385", ]
# peaking.gene <- expression.matrix["ENSG00000123405", ]
peaking.gene <- expression.matrix["ENSG00000196460", ]
peaking.gene <- melt(peaking.gene[, 1:12])
peaking.gene$tp <- hours
peaking.gene$class <- "peaking\n\n(RFX8)"

# 5. select bending gene
bending.gene <- expression.matrix["ENSG00000090104", ]
bending.gene <- melt(bending.gene[, 1:12])
bending.gene$tp <- hours
bending.gene$class <- "bending\n\n(RGS1)"

# 6. flat gene
flat.gene <- expression.matrix["ENSG00000004059", ]
flat.gene <- melt(flat.gene[, 1:12])
flat.gene$tp <- hours
flat.gene$class <- "flat\n\n(ARF5)"



genes <- rbind(upregulated.gene, peaking.gene, downregulated.gene, bending.gene, flat.gene)
genes$class <- factor(genes$class, levels = c("bending\n\n(RGS1)", 
                                              "down-\nregulated\n\n(IGLL1)", 
                                              "peaking\n\n(RFX8)", 
                                              "up-\nregulated\n\n(C1QA)", 
                                              "flat\n\n(ARF5)"))

scales_y <- list(
  `bending\n\n(RGS1)` = scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 2)),
  `down-\nregulated\n\n(IGLL1)` = scale_y_continuous(limits = c(0, 11), breaks = seq(0, 11, 2)),
  `peaking\n\n(RFX8)` = scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 2)),
  `up-\nregulated\n\n(C1QA)` = scale_y_continuous(limits = c(0, 11), breaks = seq(0, 11, 2)),
  `flat\n\n(ARF5)` = scale_y_continuous(limits = c(5, 10), breaks = seq(5, 10, 2))
)


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1h.pdf", height = 10, width=5.5)
ggplot(genes, aes(x=variable, y=value, group=1, colour=class, fill=class)) +
  geom_line(size=2) +
  facet_grid_sc(vars(class), scales = list(y = scales_y)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  theme_bw() +
  theme(axis.text.x = element_text(size=15, angle = 30, vjust = .5),
        axis.text = element_text(size=15),
        axis.title = element_text(size = 15),
        strip.text.y = element_text(size = 15, angle=0),
        strip.background.y = element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
  ) +
  guides(fill=F, color=F) +
  ylab("log2(TPM + 1)") +
  xlab("time (hours)") +
  scale_x_discrete(labels = hours)
dev.off()

