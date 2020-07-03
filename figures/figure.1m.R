.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")



# palette <- c("H3K9ac" = "#e7298a",
#              "H3K27ac" = "#8e0152",
#              "H3K4me3" = "#f1b6da",
#              "H3K27me3" = "#253494",
#              "H3K9me3" = "#41b6c4",
#              "H3K36me3" = "#7fbc41",
#              "H4K20me1" = "#276419",
#              "H3K4me1" = "#ffbf00",
#              "H3K4me2" = "#a67c00",
#              "expression" = "white")


palette <- c("H3K9ac" = "#af4d85",
             "H3K27ac" = "#630039",
             "H3K4me3" = "#d199b9",
             "H3K27me3" = "#1d2976",
             "H3K9me3" = "#a7add4",
             "H3K36me3" = "#7fbc41",
             "H4K20me1" = "#4c7027",
             "H3K4me1" = "#e5ab00",
             "H3K4me2" = "#a67c00",
             "expression" = "white")


palette2 <- c("H3K9ac" = "#af4d85",
             "H3K27ac" = "#630039",
             "H3K4me3" = "#d199b9",
             "H3K27me3" = "#1d2976",
             "H3K9me3" = "#a7add4",
             "H3K36me3" = "#7fbc41",
             "H4K20me1" = "#4c7027",
             "H3K4me1" = "#e5ab00",
             "H3K4me2" = "#a67c00",
             "expression" = "gray")




#************
# LIBRARIES *
#************

library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(cowplot)
library(reshape2)


setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/")

#********
# BEGIN *
#********

# 1. read dataframes

expression.matrix <- read.table("expression/QN.merged/selected.genes.rep.2.3.after.QN.merged.tsv",
                                h=T, sep="\t")
H3K27ac.matrix <- read.table("H3K27ac/QN.merged/H3K27ac.matrix.after.QN.merged.tsv",
                                 h=T, sep="\t")
H3K9ac.matrix <- read.table("H3K9ac/QN.merged/H3K9ac.matrix.after.QN.merged.tsv",
                            h=T, sep="\t")
H4K20me1.matrix <- read.table("H4K20me1/QN.merged/H4K20me1.matrix.after.QN.merged.tsv",
                              h=T, sep="\t")
H3K4me3.matrix <- read.table("H3K4me3/QN.merged/H3K4me3.matrix.after.QN.merged.tsv",
                             h=T, sep="\t")
H3K4me1.matrix <- read.table("H3K4me1/QN.merged/H3K4me1.matrix.after.QN.merged.tsv",
                             h=T, sep="\t")
H3K36me3.matrix <- read.table("H3K36me3/QN.merged/H3K36me3.matrix.after.QN.merged.tsv",
                              h=T, sep="\t")
H3K4me2.matrix <- read.table("H3K4me2/QN.merged/H3K4me2.matrix.after.QN.merged.tsv",
                             h=T, sep="\t")
H3K9me3.matrix <- read.table("H3K9me3/QN.merged/H3K9me3.matrix.after.QN.merged.tsv",
                             h=T, sep="\t")
H3K27me3.matrix <- read.table("H3K27me3/QN.merged/H3K27me3.matrix.after.QN.merged.tsv",
                              h=T, sep="\t")

# 2. check order of rownames
expression.matrix <- expression.matrix[rownames(H3K27ac.matrix), ]

stopifnot(identical(rownames(expression.matrix),
                    rownames(H3K27ac.matrix)))
stopifnot(identical(rownames(expression.matrix),
                    rownames(H3K9ac.matrix)))
stopifnot(identical(rownames(expression.matrix),
                    rownames(H4K20me1.matrix)))
stopifnot(identical(rownames(expression.matrix),
                    rownames(H3K4me3.matrix)))
stopifnot(identical(rownames(expression.matrix),
                    rownames(H3K4me1.matrix)))
stopifnot(identical(rownames(expression.matrix),
                    rownames(H3K36me3.matrix)))
stopifnot(identical(rownames(expression.matrix),
                    rownames(H3K4me2.matrix)))
stopifnot(identical(rownames(expression.matrix),
                    rownames(H3K9me3.matrix)))
stopifnot(identical(rownames(expression.matrix),
                    rownames(H3K27me3.matrix)))

# 3. get transpose & scaled matrices
expression.matrix <- as.data.frame(t(expression.matrix))
expression.matrix.scaled <- as.data.frame(scale(expression.matrix))
sum(apply(expression.matrix.scaled, 1, is.na))
expression.matrix.scaled <- as.data.frame(t(na.omit(t(expression.matrix.scaled))))

H3K27ac.matrix <- as.data.frame(t(H3K27ac.matrix))
H3K27ac.matrix.scaled <- as.data.frame(scale(H3K27ac.matrix))
sum(apply(H3K27ac.matrix.scaled, 1, is.na))
H3K27ac.matrix.scaled <- as.data.frame(t(na.omit(t(H3K27ac.matrix.scaled))))

H3K9ac.matrix <- as.data.frame(t(H3K9ac.matrix))
H3K9ac.matrix.scaled <- as.data.frame(scale(H3K9ac.matrix))
sum(apply(H3K9ac.matrix.scaled,1, is.na))
H3K9ac.matrix.scaled <- as.data.frame(t(na.omit(t(H3K9ac.matrix.scaled))))

H4K20me1.matrix <- as.data.frame(t(H4K20me1.matrix))
H4K20me1.matrix.scaled <- as.data.frame(scale(H4K20me1.matrix))
sum(apply(H4K20me1.matrix.scaled,1, is.na))
H4K20me1.matrix.scaled <- as.data.frame(t(na.omit(t(H4K20me1.matrix.scaled))))

H3K4me3.matrix <- as.data.frame(t(H3K4me3.matrix))
H3K4me3.matrix.scaled <- as.data.frame(scale(H3K4me3.matrix))
sum(apply(H3K4me3.matrix.scaled,1, is.na))
H3K4me3.matrix.scaled <- as.data.frame(t(na.omit(t(H3K4me3.matrix.scaled))))

H3K4me1.matrix <- as.data.frame(t(H3K4me1.matrix))
H3K4me1.matrix.scaled <- as.data.frame(scale(H3K4me1.matrix))
sum(apply(H3K4me1.matrix.scaled,1, is.na))
H3K4me1.matrix.scaled <- as.data.frame(t(na.omit(t(H3K4me1.matrix.scaled))))

H3K36me3.matrix <- as.data.frame(t(H3K36me3.matrix))
H3K36me3.matrix.scaled <- as.data.frame(scale(H3K36me3.matrix))
sum(apply(H3K36me3.matrix.scaled,1, is.na))
H3K36me3.matrix.scaled <- as.data.frame(t(na.omit(t(H3K36me3.matrix.scaled))))

H3K4me2.matrix <- as.data.frame(t(H3K4me2.matrix))
H3K4me2.matrix.scaled <- as.data.frame(scale(H3K4me2.matrix))
sum(apply(H3K4me2.matrix.scaled,1, is.na))
H3K4me2.matrix.scaled <- as.data.frame(t(na.omit(t(H3K4me2.matrix.scaled))))

H3K9me3.matrix <- as.data.frame(t(H3K9me3.matrix))
H3K9me3.matrix.scaled <- as.data.frame(scale(H3K9me3.matrix))
sum(apply(H3K9me3.matrix.scaled,1, is.na))
H3K9me3.matrix.scaled <- as.data.frame(t(na.omit(t(H3K9me3.matrix.scaled))))

H3K27me3.matrix <- as.data.frame(t(H3K27me3.matrix))
H3K27me3.matrix.scaled <- as.data.frame(scale(H3K27me3.matrix))
sum(apply(H3K27me3.matrix.scaled,1, is.na))
H3K27me3.matrix.scaled <- as.data.frame(t(na.omit(t(H3K27me3.matrix.scaled))))







###----
### PCA scaling each matrix separately
###----



# 9. PCA 

my.genes <- Reduce(intersect, list(colnames(expression.matrix.scaled),
                                   colnames(H3K27ac.matrix.scaled),
                                   colnames(H3K9me3.matrix.scaled),
                                   colnames(H4K20me1.matrix.scaled),
                                   colnames(H3K4me3.matrix.scaled),
                                   colnames(H3K4me1.matrix.scaled),
                                   colnames(H3K36me3.matrix.scaled),
                                   colnames(H3K4me2.matrix.scaled),
                                   colnames(H3K9me3.matrix.scaled),
                                   colnames(H3K27me3.matrix.scaled)))


pca.all.marks.2 <- prcomp(rbind(expression.matrix.scaled[, my.genes],
                                H3K27ac.matrix.scaled[, my.genes],
                                H3K9ac.matrix.scaled[, my.genes],
                                H4K20me1.matrix.scaled[, my.genes],
                                H3K4me3.matrix.scaled[, my.genes],
                                H3K4me1.matrix.scaled[, my.genes],
                                H3K36me3.matrix.scaled[, my.genes],
                                H3K4me2.matrix.scaled[, my.genes],
                                H3K9me3.matrix.scaled[, my.genes],
                                H3K27me3.matrix.scaled[, my.genes]),
                          center = F,
                          scale = F)

summary(pca.all.marks.2)

df.pca.all.marks.2 <- as.data.frame(pca.all.marks.2$x[, 1:4])
df.pca.all.marks.2$tp_type <- c(paste0(rownames(df.pca.all.marks.2)[1:12], "_expression"),
                                paste0(rownames(df.pca.all.marks.2)[1:12], "_H3K27ac"),
                                paste0(rownames(df.pca.all.marks.2)[1:12], "_H3K9ac"),
                                paste0(rownames(df.pca.all.marks.2)[1:12], "_H4K20me1"),
                                paste0(rownames(df.pca.all.marks.2)[1:12], "_H3K4me3"),
                                paste0(rownames(df.pca.all.marks.2)[1:12], "_H3K4me1"),
                                paste0(rownames(df.pca.all.marks.2)[1:12], "_H3K36me3"),
                                paste0(rownames(df.pca.all.marks.2)[1:12], "_H3K4me2"),
                                paste0(rownames(df.pca.all.marks.2)[1:12], "_H3K9me3"),
                                paste0(rownames(df.pca.all.marks.2)[1:12], "_H3K27me3"))
df.pca.all.marks.2$tp <- rep(rownames(df.pca.all.marks.2)[1:12], 10)
df.pca.all.marks.2$type <- rep(c("expression", "H3K27ac",
                                 "H3K9ac", "H4K20me1", "H3K4me3",
                                 "H3K4me1", "H3K36me3",
                                 "H3K4me2", "H3K9me3",
                                 "H3K27me3"), each=12)

df.pca.all.marks.2$type <- factor(df.pca.all.marks.2$type, levels = c("expression",
                                                                      "H3K27ac",
                                                                      "H3K9ac",
                                                                      "H4K20me1",
                                                                      "H3K4me3",
                                                                      "H3K4me1",
                                                                      "H3K36me3",
                                                                      "H3K4me2",
                                                                      "H3K9me3",
                                                                      "H3K27me3"))

df.pca.all.marks.2$tp_numeric <- rep(c(0,3,6,9,12,18,24,36,48,72,120,168), 10)
df.pca.all.marks.2$tp_cat <- rep(1:12, 10)



# 10. PC1 vs PC2


p <- ggplot(df.pca.all.marks.2,
       aes(x=PC1, y=PC2, fill = type, label = tp_numeric, colour=type)) +
  theme_bw() +
  facet_wrap(~type, nrow=2) +
  geom_density_2d(colour="white") + 
  theme(axis.title = element_text(size =13),
        axis.text = element_text(size = 13),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        # aspect.ratio = 1,
        strip.background.x = element_blank(),
        strip.text.x = element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.position = "bottom",
        plot.title = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  geom_point(data = df.pca.all.marks.2,
             aes(size=tp_numeric), shape=21, alpha = .8, colour="black") +
  geom_line(data = df.pca.all.marks.2, alpha=.3, size=10) +
  #geom_text_repel(data = df.pca.all.marks.2) +
  scale_fill_manual(values = palette) +
  scale_colour_manual(values = palette2) +
  xlab("PC1 (13.8 %)") +
  ylab("PC2 (5.1 %)") +
  guides(fill=F,
         colour=F,
         size = guide_legend(title = "hours"))

my.legend <- get_legend(p)
p <- p + guides(size=F)

main.p <- plot_grid(plotlist = list(p, my.legend),
                    nrow = 2, rel_heights = c(1, 0.1))

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1m.pdf", width = 12, height = 6.5)
print(main.p)
dev.off()
