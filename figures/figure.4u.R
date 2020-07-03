.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(reshape2)
library(ggalluvial)
library(viridis)
library(ggpubr)


#**********
# PALETTE *
#**********

palette <- c("exp_ant" = "black",
             "mark_ant" = "#95A10D",
             "no_diff" = "#bbbbbb")




#************
# FUNCTIONS *
#************

my.function1 <- function(my.df, my.vector, my.column) {
  
  for ( i in 1:nrow(my.df) ){
    
    if (my.df[i, my.column] < 0 ) {
      
      my.vector <- c(my.vector, "exp_ant")
      
    } else if (my.df[i, my.column] > 0 ) {
      
      my.vector <- c(my.vector, "mark_ant")
      
    } else {
      
      my.vector <- c(my.vector, "no_diff")
      
    }
    
  } 
  
  my.df2 <- data.frame(group=my.vector,
                       gene_id=rownames(my.df))
  
  return(my.df2)
  
}



my.function2 <- function(mark) {
  
  # 1. import expression df for time-points
  # where 25%, 50%, 75%, 100% of up-regulation occurs
  expression.m <- read.table(paste0(mark, "/QN.merged/ant.del.analysis/increases/expression.25.50.75.100.tsv"),
                             h=T, sep="\t")
  expression.m$gene_id <- rownames(expression.m)
  
  
  # 2. import mark df for time-points where
  # 25%, 50%, 75%, 100% of up-regulation occurs
  mark.m <- read.table(paste0(mark, "/QN.merged/ant.del.analysis/increases/", mark, ".25.50.75.100.tsv"),
                       h=T, sep="\t")
  mark.m$gene_id <- rownames(mark.m)
  
  
  # 3. import list of positively_correlated & up-regulated genes 
  upreg.pc <- read.table(paste0(mark, "/QN.merged/ant.del.analysis/increases/upregulation.positively_correlated.genes.txt"))
  
  
  # 4. merge the expression and mark dataframes
  merged.m <- merge(expression.m, mark.m, by = "gene_id")
  
  # 5. keep only positively_correlated & up-regulated genes 
  merged.m <- merged.m[merged.m$gene_id %in% upreg.pc$V1, ]
  rownames(merged.m) <- merged.m$gene_id
  merged.m$gene_id <- NULL
  colnames(merged.m) <- paste(rep(c("expression", "mark"), each=4),
                              rep(c(25, 50, 75, 100), 2), sep="_")
  
  
  # 6. compute, for 25% of up-regulation,
  # the difference between the time-points of expression and the mark
  merged.m$perc_25 <- merged.m$expression_25 - merged.m$mark_25
  merged.m2 <- merged.m[, "perc_25", drop=F]
  
  
  # 7. assign to each gene and for 25% of up-regulation
  # one of these categories: exp_ant, mark_ant, no_diff
  
  merged.m3.perc_25 <- c()
  merged.m3.perc_25 <- my.function1(my.df=merged.m2, my.column="perc_25", my.vector=merged.m3.perc_25)
  
  
  # 8. read original expression matrix (log2 TPM)
  my.expression <- read.table(paste0(mark, "/QN.merged/", mark, ".matrix.tsv"),
                              h=T, sep="\t")
  
  my.expression <- my.expression[, "H000", drop=F]
  my.expression$gene_id <- rownames(my.expression)
  
  
  # 9. retrieve expression at 0h for genes belonging to exp_ant, mark_ant, no_diff
  merged.m3.perc_25 <- merge(merged.m3.perc_25,
                             my.expression,
                             by = "gene_id")
  print(mark)
  print(nrow(merged.m3.perc_25))
  
  # 10. prepare df for alluvial plot
  merged.m3.perc_25$group <- factor(merged.m3.perc_25$group, 
                                    levels=c("exp_ant", "no_diff", "mark_ant"))
  merged.m3.perc_25 <- merged.m3.perc_25[merged.m3.perc_25$group != "no_diff", ]
  my.max <- max(merged.m3.perc_25$H000) + 0.25
  
  # 11. make boxplot 
  p <- ggplot(merged.m3.perc_25,
              aes(x = group, y=H000, fill=group)) +
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
    ylab("") +
    stat_compare_means(comparisons = list(c("exp_ant", "mark_ant")),
                       method = "wilcox.test",
                       size = 4.5,
                       label.y = my.max - 0.1) +
    scale_fill_manual(values = palette) +
    labs(title=mark) +
    ylim(0, my.max)
  
  return(p)
  
  
}


#********
# BEGIN *
#********

lop <- list()

marks <- c("H3K27ac", "H3K4me1", "H3K4me2", "H3K9ac", "H3K4me3", "H3K36me3", "H4K20me1")

for ( i in 1:7 ){
  
  lop[[i]] <- my.function2(mark=marks[i])
  
}


my.legend <- get_legend(lop[[1]])

lop <- lapply(lop, function(x){x <- x+guides(fill=F, color=F)})
lop[[8]] <- my.legend
lop[[1]] <- lop[[1]] + ylab("mark level at 0h")

dev.off()

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.4u.pdf", width=15, height=3.5)
plot_grid(plotlist = lop, nrow=1, ncol=7)
dev.off()
