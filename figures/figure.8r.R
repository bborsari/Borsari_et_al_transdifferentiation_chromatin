.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(reshape2)
library(ggalluvial)
library(viridis)
library(ggpubr)

#************
# FUNCTIONS *
#************

function1 <- function(mark) {
  
  df <- read.table(paste0(mark, "/QN.merged/ant.del.analysis/increases/", 
                          "expression.", mark, ".dynamics.tsv"),
                   h=T, sep="\t")
  
  df <- df[, c("group_25", "group_50", "group_75", "group_100")]
  
  df$gene_id <- rownames(df)
  
  df.melt <- melt(df, id.vars = "gene_id")
  df.melt$value <- factor(df.melt$value, levels = c("delayed", "concomitant", "anticipated"))
  
  df.melt <- df.melt[df.melt$variable == "group_25" & df.melt$value == "anticipated", ]
  
  df.melt <- merge(df.melt, 
                   expression.m[, c("gene_id", "H000")],
                   by = "gene_id")
  
  df.melt$mark <- mark
  
  return(df.melt)
  
}

#********
# BEGIN *
#********

# 1. set working directory
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


# 2. define palette
palette <- c("H3K9ac" = "#e7298a",
             "H3K27ac" = "#8e0152",
             "H3K4me3" = "#f1b6da",
             "H3K27me3" = "#253494",
             "H3K9me3" = "#41b6c4",
             "H3K36me3" = "#7fbc41",
             "H4K20me1" = "#276419",
             "H3K4me1" = "#ffbf00",
             "H3K4me2" = "#a67c00")


# 3. the marks we're analyzing
marks <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3", "H3K36me3", "H4K20me1")


# 4. read expression matrix
expression.m <- read.table("H3K4me3/QN.merged/expression.matrix.tsv", h=T, sep="\t",
                           stringsAsFactors = F)
expression.m$gene_id <- rownames(expression.m)


# 5. build a unique df for all marks
x <- data.frame(stringsAsFactors = F)
for ( i in 1:7 ) {
  
  x <- rbind(x, function1(mark = marks[i]))
  
}

x$mark <- factor(x$mark, levels = marks)


# 6. make plot
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8r.pdf", 
    width=7, height=3.5)
ggplot(x, aes(x = mark, y=H000, fill=mark)) +
  geom_violin(alpha=.4, colour="white") +
  geom_boxplot(width=0.25, alpha = .4, outlier.shape = NA) +
  guides(fill=F) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=11),
        plot.title = element_text(size = 17, hjust = .5),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        axis.text.y = element_text(size =11, angle=90, hjust=.5),
        axis.text.x = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        legend.title = element_blank()) +
  ylab("expression at 0h - log2 (TPM+1)") +
  scale_fill_manual(values = palette) +
  ylim(0, 12)
dev.off()
