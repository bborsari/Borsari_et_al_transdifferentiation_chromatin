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


#**********
# PALETTE *
#**********



# palette <- c("mark_ant" = "#a1d99b",
#              "no_diff" = "#41ab5d",
#              "exp_ant" = "#00441b")

palette <- c("mark_ant" = "#bdbdbd",
             "no_diff" = "#737373",
             "exp_ant" = "#000000")



#************
# FUNCTIONS *
#************

my.function1 <- function(f1.my.df, f1.my.vector, f1.my.column) {
  
  for ( i in 1:nrow(f1.my.df) ){
    
    if (f1.my.df[i, f1.my.column] < 0 ) {
      
      f1.my.vector <- c(f1.my.vector, "exp_ant")
      
    } else if (f1.my.df[i, f1.my.column] > 0 ) {
      
      f1.my.vector <- c(f1.my.vector, "mark_ant")
      
    } else {
      
      f1.my.vector <- c(f1.my.vector, "no_diff")
      
    }
    
  } 
  
  return(f1.my.vector)
  
}


my.function2 <- function(f2.mark) {
  
  # 1. import expression df for time-points
  # where 25%, 50%, 75%, 100% of up-regulation occurs
  f2.expression.m <- read.table(paste0(f2.mark, "/QN.merged/ant.del.analysis/increases/expression.25.50.75.100.tsv"),
                                h=T, sep="\t")
  f2.expression.m$gene_id <- rownames(f2.expression.m)
  
  
  # 2. import mark df for time-points where
  # 25%, 50%, 75%, 100% of up-regulation occurs
  f2.mark.m <- read.table(paste0(f2.mark, "/QN.merged/ant.del.analysis/increases/", 
                                 f2.mark, ".25.50.75.100.tsv"),
                          h=T, sep="\t")
  f2.mark.m$gene_id <- rownames(f2.mark.m)
  
  
  # 3. import list of positively_correlated & up-regulated genes 
  f2.upreg.pc <- read.table(paste0(f2.mark, 
                                   "/QN.merged/ant.del.analysis/increases/upregulation.positively_correlated.genes.txt"))
  
  
  # 4. merge the expression and mark dataframes
  f2.merged.m <- merge(f2.expression.m, f2.mark.m, by = "gene_id")
  
  # 5. keep only positively_correlated & up-regulated genes 
  f2.merged.m <- f2.merged.m[f2.merged.m$gene_id %in% f2.upreg.pc$V1, ]
  rownames(f2.merged.m) <- f2.merged.m$gene_id
  f2.merged.m$gene_id <- NULL
  colnames(f2.merged.m) <- paste(rep(c("expression", "mark"), each=4),
                                 rep(c(25, 50, 75, 100), 2), sep="_")
  
  
  # 6. compute, for each degree of up-regulation (25%, 50%, 75%, 100%),
  # the difference between the time-points of expression and the mark
  f2.merged.m$perc_25 <- f2.merged.m$expression_25 - f2.merged.m$mark_25
  f2.merged.m$perc_50 <- f2.merged.m$expression_50 - f2.merged.m$mark_50
  f2.merged.m$perc_75 <- f2.merged.m$expression_75 - f2.merged.m$mark_75
  f2.merged.m$perc_100 <- f2.merged.m$expression_100 - f2.merged.m$mark_100
  f2.merged.m2 <- f2.merged.m[, c("perc_25", "perc_50", "perc_75", "perc_100")]
  
  
  # 7. assign to each gene and to each degree of up-regulation
  # one of these categories: exp_ant, mark_ant, no_diff
  
  f2.merged.m3.perc_25 <- c()
  f2.merged.m3.perc_50 <- c()
  f2.merged.m3.perc_75 <- c()
  f2.merged.m3.perc_100 <- c()
  
  f2.merged.m3.perc_25 <- my.function1(f1.my.df=f2.merged.m2, 
                                       f1.my.column="perc_25", 
                                       f1.my.vector=f2.merged.m3.perc_25)
  f2.merged.m3.perc_50 <- my.function1(f1.my.df=f2.merged.m2, 
                                       f1.my.column="perc_50", 
                                       f1.my.vector=f2.merged.m3.perc_50)
  f2.merged.m3.perc_75 <- my.function1(f1.my.df=f2.merged.m2, 
                                       f1.my.column="perc_75", 
                                       f1.my.vector=f2.merged.m3.perc_75)
  f2.merged.m3.perc_100 <- my.function1(f1.my.df=f2.merged.m2, 
                                        f1.my.column="perc_100", 
                                        f1.my.vector=f2.merged.m3.perc_100)
  
  print(f2.mark)
  print(table(f2.merged.m3.perc_25)/length(f2.merged.m3.perc_25))
  
  # 8. prepare df for alluvial plot
  f2.merged.m3 <- data.frame(f2.merged.m3.perc_25,
                             f2.merged.m3.perc_50,
                             f2.merged.m3.perc_75,
                             f2.merged.m3.perc_100)
  f2.merged.m3$gene_id <- rownames(f2.merged.m2)
  f2.merged.m3 <- melt(f2.merged.m3, id.vars = "gene_id")
  f2.merged.m3$value <- factor(f2.merged.m3$value, levels=c("exp_ant", "no_diff", "mark_ant"))
  
  
  # 9. make alluvial plot  
  f2.p <- ggplot(f2.merged.m3,
                 aes(x = variable, stratum = value, alluvium = gene_id, label = value)) +
    # geom_flow(stat = "alluvium", lode.guidance = "forward", aes(fill=value), alpha=.6) +
    geom_flow(aes(fill = value), alpha=0.8) +
    geom_stratum(aes(fill = value), alpha=0.8, color="white") +
    theme_bw() + 
    theme(panel.border = element_rect(color="black"),
          plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_blank(),
          axis.title.y = element_text(size=11),
          axis.text.y = element_text(size=11, angle=90, hjust=.5),
          axis.text.x = element_text(size=12),
          axis.title.x = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size=15),
          legend.title = element_blank(),
          plot.title = element_text(size=17, hjust=.5)) +
    ylab("") +
    scale_x_discrete(labels = c("25%", "50%", "75%", "100%")) +
    scale_color_manual(values=palette) +
    scale_fill_manual(values=palette) +
    labs(title = f2.mark)
  
  return(f2.p)
  
  
}


#********
# BEGIN *
#********

lop <- list()

marks <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3", "H3K36me3", "H4K20me1")

for ( i in 1:7 ){
  
  lop[[i]] <- my.function2(f2.mark=marks[i])
  
}


my.legend <- get_legend(lop[[1]])

lop <- lapply(lop, function(x){x <- x+guides(fill=F, color=F)})
lop[[8]] <- my.legend
lop[[1]] <- lop[[1]] + ylab("number of genes")


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.4p.pdf", width=15, height=3)
p1 <- plot_grid(plotlist = lop[1:7], nrow=1, ncol=7)
plot_grid(plotlist = list(p1, lop[[8]]), nrow=2, rel_heights = c(1, 0.1))
dev.off()
