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

hours <- c(0, 3, 6, 9, 12, 18, 24, 36, 48, 72, 120, 168)
names(hours) <- c("H000", "H003", "H006", "H009", "H012",
                  "H018", "H024", "H036", "H048", "H072",
                  "H120", "H168")

my.function <- function(mark, type="pearson") {
  
  mark.matrix <- read.table(paste0(mark, "/QN.merged/", mark, ".matrix.after.QN.merged.tsv"),
                            h=T, sep="\t")
  
  stopifnot(identical(rownames(expression.matrix),
                      rownames(mark.matrix)))
  
  mark.cc <- diag(cor(expression.matrix, mark.matrix, method=type))
  mark.cc.df <- data.frame(cc=mark.cc)
  mark.cc.df$mark <- mark
  mark.cc.df$tp <- rownames(mark.cc.df)
  
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

all.marks.df.pearson$tp <- hours[all.marks.df.pearson$tp]


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.10e.pearson.pdf",
    width = 10,
    height = 3)
ggplot(all.marks.df.pearson, 
       aes(x=mark, y=cc, fill=mark, size=tp)) +
  geom_jitter(shape=21) +
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
  ) +
  ylim(-1, 1)
dev.off()


