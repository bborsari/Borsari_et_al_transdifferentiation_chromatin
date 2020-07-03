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
  df.melt$value <- factor(df.melt$value, levels = c("anticipated", "concomitant", "delayed"))
  
  df.melt <- df.melt[df.melt$variable == "group_25" & 
                       df.melt$value != "concomitant", ]
  
  df.melt <- merge(df.melt, 
                   expression.m[, c("gene_id", "H168")],
                   by = "gene_id")
  
  df.melt$mark <- mark
  
  p <- ggplot(df.melt, aes(x = value, y=H168, fill=value)) +
    geom_violin(alpha=.6, colour="white") +
    geom_boxplot(width=0.25, alpha = .4, outlier.shape = NA) +
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
          legend.title = element_blank(),
          legend.text = element_text(size=15),
          legend.position = "bottom") +
    ylab("expression at 168h - log2 (TPM+1)") +
    scale_fill_manual(values = palette) +
    ylim(0, 16) +
    stat_compare_means(method = "wilcox.test",
                       size = 4.5,
                       label.y = 15,
                       label.x = 1.2,
                       label = "p.format") +
    labs(title = mark)
  
  
}

#********
# BEGIN *
#********

# 1. set working directory
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


# 2. define palette
palette <- c("anticipated" = "#bdbdbd",
             "delayed" = "#403734")

# 3. the marks we're analyzing
marks <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3", "H3K36me3", "H4K20me1")


# 4. read expression matrix
expression.m <- read.table("H3K4me3/QN.merged/expression.matrix.tsv", h=T, sep="\t",
                           stringsAsFactors = F)
expression.m$gene_id <- rownames(expression.m)


# 5. make plots
lop <- list()
for ( i in 1:7 ) {
  
  lop[[i]] <- function1(mark = marks[i])
  if ( i != 1 ){
    
    lop[[i]] <- lop[[i]] + ylab("")
    
  }
  
}


leg <- get_legend(lop[[1]])
lop <- lapply(lop, function(x){x <- x + 
  guides(fill=F) +
  scale_x_discrete(labels = c("a", "d"))})



# 6. make plot
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.9x.pdf", 
    width=15, height=3.5)
main <- plot_grid(plotlist = lop, nrow=1, align="h")
plot_grid(main, leg, nrow=2, rel_heights = c(1, 0.1))
dev.off()
