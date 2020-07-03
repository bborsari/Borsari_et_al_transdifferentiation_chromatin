.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/")
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



#************
# LIBRARIES *
#************

library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(cowplot)
library(reshape2)
library(dplyr)
library(tidyr)
library(ggcorrplot)
library(ggsignif)
library(ggpubr)



#************
# FUNCTIONS *
#************

my.function <- function(mark, type="pearson") {
  
  mark.matrix <- read.table(paste0(mark, "/QN.merged/", mark, ".matrix.after.QN.merged.tsv"),
                                   h=T, sep="\t")
  
  stopifnot(identical(rownames(expression.matrix),
                      rownames(mark.matrix)))
  
  mark.cc <- diag(cor(t(expression.matrix), t(mark.matrix), method=type))
  mark.cc.df <- data.frame(cc=mark.cc)
  rownames(mark.cc.df) <- rownames(expression.matrix)
  mark.cc.df$gene_id <- rownames(mark.cc.df)
  mark.cc.df$mark <- mark
  
  return(mark.cc.df)


}

#********
# BEGIN *
#********

# 1. read expression matrix
expression.matrix <- read.table("expression/QN.merged/selected.genes.rep.2.3.after.QN.merged.tsv",
                                h=T, sep="\t")
# 2. sort rownames of expression matrix
new.order <- sort(rownames(expression.matrix))
expression.matrix <- expression.matrix[new.order, ]


# 3. compute cc for all marks 
all.marks.df.pearson <- data.frame()

for (my.mark in c("H3K27ac", "H3K9ac", "H4K20me1",
               "H3K36me3", "H3K4me3", "H3K4me1",
               "H3K4me2", "H3K9me3", "H3K27me3")) {
  
  all.marks.df.pearson <- rbind(all.marks.df.pearson,
                                my.function(mark=my.mark))
  
}


# 4. compute mean, median
means.pearson <- c()
medians.pearson <- c()

for (my.mark in c("H3K27ac", "H3K9ac", "H4K20me1",
                  "H3K36me3", "H3K4me3", "H3K4me1",
                  "H3K4me2", "H3K9me3", "H3K27me3")) {
  
  means.pearson <- c(means.pearson, 
                     mean(all.marks.df.pearson[all.marks.df.pearson$mark == my.mark, "cc"], 
                          na.rm = T))
  medians.pearson <- c(medians.pearson, 
                       median(all.marks.df.pearson[all.marks.df.pearson$mark == my.mark, "cc"], 
                              na.rm = T))
  
  
}

names(means.pearson) <- c("H3K27ac", "H3K9ac", "H4K20me1",
                  "H3K36me3", "H3K4me3", "H3K4me1",
                  "H3K4me2", "H3K9me3", "H3K27me3")
names(medians.pearson) <- c("H3K27ac", "H3K9ac", "H4K20me1",
                    "H3K36me3", "H3K4me3", "H3K4me1",
                    "H3K4me2", "H3K9me3", "H3K27me3")

names(means.pearson[order(means.pearson, decreasing = T)])
names(medians.pearson[order(medians.pearson, decreasing = T)])

all.marks.df.pearson$mark <- as.factor(all.marks.df.pearson$mark)
all.marks.df.pearson$mark <- factor(all.marks.df.pearson$mark, 
                                    levels = c("H3K27ac",
                                               "H3K9ac",
                                               "H4K20me1",
                                               "H3K4me3",
                                               "H3K4me1",
                                               "H3K36me3",
                                               "H3K4me2",
                                               "H3K9me3",
                                               "H3K27me3"))

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.2c.pearson.pdf",
    width = 10,
    height = 3)
ggplot(all.marks.df.pearson, 
       aes(x=mark, y=cc, fill=mark)) +
  geom_violin(alpha=.5, width=.7, colour="white") +
  geom_boxplot(width=0.2, alpha = .6, outlier.shape = NA) +
  guides(fill=F) +
  ylab("Pearson's R") +
  scale_fill_manual(values = palette) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text = element_text(size=13),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")
  ) 
dev.off()


#*************************

# 3. compute cc for all marks 
all.marks.df.spearman <- data.frame()

for (my.mark in c("H3K27ac", "H3K9ac", "H4K20me1",
                  "H3K36me3", "H3K4me3", "H3K4me1",
                  "H3K4me2", "H3K9me3", "H3K27me3")) {
  
  all.marks.df.spearman <- rbind(all.marks.df.spearman,
                                my.function(mark=my.mark, type="spearman"))
  
}


# 4. compute mean, median
means.spearman <- c()
medians.spearman <- c()

for (my.mark in c("H3K27ac", "H3K9ac", "H4K20me1",
                  "H3K36me3", "H3K4me3", "H3K4me1",
                  "H3K4me2", "H3K9me3", "H3K27me3")) {
  
  means.spearman <- c(means.spearman, 
                     mean(all.marks.df.spearman[all.marks.df.spearman$mark == my.mark, "cc"], 
                          na.rm = T))
  medians.spearman <- c(medians.spearman, 
                       median(all.marks.df.spearman[all.marks.df.spearman$mark == my.mark, "cc"], 
                              na.rm = T))
  
  
}

names(means.spearman) <- c("H3K27ac", "H3K9ac", "H4K20me1",
                          "H3K36me3", "H3K4me3", "H3K4me1",
                          "H3K4me2", "H3K9me3", "H3K27me3")
names(medians.spearman) <- c("H3K27ac", "H3K9ac", "H4K20me1",
                            "H3K36me3", "H3K4me3", "H3K4me1",
                            "H3K4me2", "H3K9me3", "H3K27me3")

names(means.spearman[order(means.spearman, decreasing = T)])
names(medians.spearman[order(medians.spearman, decreasing = T)])

all.marks.df.spearman$mark <- as.factor(all.marks.df.spearman$mark)
all.marks.df.spearman$mark <- factor(all.marks.df.spearman$mark, 
                                    levels = c("H3K27ac",
                                               "H3K9ac",
                                               "H4K20me1",
                                               "H3K4me3",
                                               "H3K4me1",
                                               "H3K36me3",
                                               "H3K4me2",
                                               "H3K9me3",
                                               "H3K27me3"))

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.2c.spearman.pdf",
    width = 10,
    height = 3)
ggplot(all.marks.df.spearman, 
       aes(x=mark, y=cc, fill=mark)) +
  geom_violin(alpha=.65) +
  geom_boxplot(width=0.3, alpha = .85, outlier.shape = NA) +
  guides(fill=F) +
  ylab("spearman's R") +
  scale_fill_manual(values = palette) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text = element_text(size=13),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")
  ) 
dev.off()


