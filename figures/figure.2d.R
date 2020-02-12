.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")



palette <- c("down-\nregulated" = "#810f7c", 
             "bending" = "#737373",
             "up-\nregulated" = "#f16913",
             "peaking" = "#c7e9b4",
             "flat" = "#1c9099")


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
library(ggridges)
library(lemon)


#************
# FUNCTIONS *
#************

my.function <- function(mark, type="pearson") {
  
  mark.matrix <- read.table(paste0(mark, "/QN.merged/", mark, ".matrix.after.QN.merged.tsv"),
                                   h=T, sep="\t")
  
  mark.matrix <- mark.matrix[rownames(mark.matrix) %in% rownames(metadata), ]
  mark.matrix <- mark.matrix[rownames(metadata), ]
  
  stopifnot(identical(rownames(expression.matrix),
                      rownames(mark.matrix)))
  
  mark.cc <- diag(cor(t(expression.matrix), t(mark.matrix), method=type))
  mark.cc.df <- data.frame(cc=mark.cc, gene_id=rownames(expression.matrix))
  mark.cc.df <- merge(mark.cc.df, metadata, by = "gene_id")
  mark.cc.df$mark <- mark
  return(mark.cc.df)
  
}


#********
# BEGIN *
#********

# 1. read dataframes
expression.matrix <- read.table("expression/QN.merged/selected.genes.rep.2.3.after.QN.merged.tsv",
                                h=T, sep="\t")
metadata <- read.table("expression/QN.merged/metadata.class2.tsv", h=T, sep="\t")
metadata$class <- gsub("regulation", "-\nregulated", metadata$class)
metadata$class2 <- gsub("regulation", "-\nregulated", metadata$class2)
metadata$class3 <- ifelse(metadata$class %in% c("down-\nregulated", "up-\nregulated", "peaking", "bending"),
                          as.character(metadata$class), paste0(metadata$class2, "#2"))
metadata$class4 <- ifelse(metadata$class %in% c("down-\nregulated", "up-\nregulated", "peaking", "bending"),
                          as.character(metadata$class), as.character(metadata$class2))
silent.genes <- read.table("expression/silent.genes.txt", h=F, sep="\t")
silent.genes$V1 <- gsub("\\..*", "", silent.genes$V1)
stable.genes <- setdiff(rownames(expression.matrix), rownames(metadata))
stable.genes <- setdiff(stable.genes, silent.genes$V1)
stable.genes.metadata <- data.frame(class = rep("flat", length(stable.genes)),
                                    time_point = rep(NA, length(stable.genes)),
                                    hc = rep(NA, length(stable.genes)),
                                    class2 = rep(NA, length(stable.genes)),
                                    class3 = rep(NA, length(stable.genes)),
                                    class4 = rep("flat", length(stable.genes)))
rownames(stable.genes.metadata) <- stable.genes
metadata <- rbind(metadata, stable.genes.metadata)
expression.matrix <- expression.matrix[rownames(expression.matrix) %in% rownames(metadata), ]
expression.matrix <- expression.matrix[rownames(metadata), ]
stopifnot(identical(rownames(expression.matrix), rownames(metadata)))
metadata$gene_id <- rownames(metadata)

# 2. distribution of Pearson's cc across time for the different marks

## H3K27ac
H3K27ac.df <- my.function(mark = "H3K27ac")

## H3K9ac
H3K9ac.df <- my.function(mark = "H3K9ac")

## H4K20me1
H4K20me1.df <- my.function(mark = "H4K20me1")

## H3K4me3
H3K4me3.df <- my.function(mark = "H3K4me3")

## H3K4me1
H3K4me1.df <- my.function(mark = "H3K4me1")

## H3K36me3
H3K36me3.df <- my.function(mark = "H3K36me3")

## H3K4me2
H3K4me2.df <- my.function(mark = "H3K4me2")

## H3K9me3
H3K9me3.df <- my.function(mark = "H3K9me3")

## H3K27me3
H3K27me3.df <- my.function(mark = "H3K27me3")


# 3. plot

## first plot
all.marks <- rbind(H3K27ac.df, H3K9ac.df, H4K20me1.df, 
                   H3K4me3.df, H3K4me1.df, H3K36me3.df, 
                   H3K4me2.df, H3K9me3.df, H3K27me3.df)

all.marks$mark <- factor(all.marks$mark,
                         levels = c("H3K27ac",
                                    "H3K9ac",
                                    "H4K20me1",
                                    "H3K4me3",
                                    "H3K4me1",
                                    "H3K36me3",
                                    "H3K4me2",
                                    "H3K9me3",
                                    "H3K27me3"))

all.marks$class4 <- factor(all.marks$class4,
                           levels = c("up-\nregulated",
                                      "peaking",
                                      "bending",
                                      "down-\nregulated",
                                      "flat"))


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.2d.pdf", width=18, height=7)
ggplot(all.marks, aes(x = cc, y = class4, fill = class4, vline_color = ..quantile.. )) + 
  geom_density_ridges(alpha = .6, quantile_lines = TRUE, quantiles = 2) +
  theme_ridges() + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 20, hjust = .5),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=16, angle = 30, vjust = .5),
        strip.text.x = element_text(size=16))  +
  scale_fill_manual(values = palette) +
  facet_grid(~mark) +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1)) +
  scale_discrete_manual("vline_color", values =c("black", "white")) +
  guides(vline_color = F) +
  xlab("Pearson's R")
dev.off()
