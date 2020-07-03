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
library(ggridges)
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

hours <- c("1" = 0, "2" = 3, "3" = 6, "4" = 9,
           "5" = 12, "6" = 18, "7" = 24,
           "8" = 36, "9" = 48, "10" = 72,
           "11" = 120, "12" = 168)

my.function1 <- function(degree, expression.mark.m, mark) {
  
  my.list <- list(c(1,2,6), c(1,3,7), c(1,4,8), c(1,5,9))
  names(my.list) <- c("perc_25", "perc_50", "perc_75", "perc_100")
  
  expression.mark.m.subset <- expression.mark.m[, my.list[[degree]]]
  colnames(expression.mark.m.subset) <- c("gene_id", "expression", "mark")
  expression.mark.m.subset.melt <- melt(expression.mark.m.subset)

  
  my.out <- expression.mark.m.subset.melt
  my.out$degree <- degree
  my.out$my.mark <- mark
  
  return(my.out)
  
}


my.function2 <- function(mark) {
  
  
  
  # 1. read expression matrix
  expression.m <- read.table(paste0(mark, "/QN.merged/ant.del.analysis/increases/expression.25.50.75.100.tsv"),
                             h=T, sep="\t")
  expression.m$gene_id <- rownames(expression.m)
  
  # 2. read mark matrix
  mark.m <- read.table(paste0(mark, "/QN.merged/ant.del.analysis/increases/", mark, ".25.50.75.100.tsv"),
                       h=T, sep="\t")
  mark.m$gene_id <- rownames(mark.m)
  
  # 3. merge mark and expression matrices
  expression.mark.m <- merge(expression.m, mark.m, by="gene_id")
  
  # 4. keep only up-regulated and positively_correlated genes
  upreg.pos.cc <- read.table(paste0(mark, "/QN.merged/ant.del.analysis/increases/upregulation.positively_correlated.genes.txt"), 
                             h=F, sep="\t", stringsAsFactors = F)
  upreg.pos.cc <- upreg.pos.cc$V1
  expression.mark.m <- expression.mark.m[expression.mark.m$gene_id %in% upreg.pos.cc, ]
  
  expression.mark.m[, 2:ncol(expression.mark.m)] <- t(apply(expression.mark.m[, 2:ncol(expression.mark.m)], 1, function(x){x <- hours[x]}))
  

  # 6. prepare df for plot
  plot.df <- rbind(my.function1(degree="perc_25", 
                                expression.mark.m = expression.mark.m,
                                mark = mark),
                   my.function1(degree="perc_50", 
                                expression.mark.m = expression.mark.m,
                                mark = mark),
                   my.function1(degree="perc_75", 
                                expression.mark.m = expression.mark.m,
                                mark = mark),
                   my.function1(degree="perc_100", 
                                expression.mark.m = expression.mark.m,
                                mark = mark))
  
  return(plot.df)
  
  
  
}


my.function3 <- function(type) {
  
  if (type == "expression") {
    
    x <- unique(my.plot.df[my.plot.df$degree == "perc_25" & my.plot.df$variable == type & my.plot.df$value == 120, "gene_id"])
    print(unique(my.plot.df[my.plot.df$gene_id %in% x & my.plot.df$degree == "perc_50" & my.plot.df$variable == type, "value"]))
    print(unique(my.plot.df[my.plot.df$gene_id %in% x & my.plot.df$degree == "perc_75" & my.plot.df$variable == type, "value"]))
    print(unique(my.plot.df[my.plot.df$gene_id %in% x & my.plot.df$degree == "perc_100" & my.plot.df$variable == type, "value"]))
    
  } else {
    
    x <- unique(my.plot.df[my.plot.df$degree == "perc_25" & 
                             my.plot.df$variable == "mark" & 
                             my.plot.df$my.mark == type & 
                             my.plot.df$value == 120, "gene_id"])
    
    print(unique(my.plot.df[my.plot.df$gene_id %in% x & 
                              my.plot.df$degree == "perc_50" & 
                              my.plot.df$variable == "mark" & 
                              my.plot.df$my.mark == type, "value"]))
    
    print(unique(my.plot.df[my.plot.df$gene_id %in% x & 
                              my.plot.df$degree == "perc_75" & 
                              my.plot.df$variable == "mark" & 
                              my.plot.df$my.mark == type, "value"]))
    
    print(unique(my.plot.df[my.plot.df$gene_id %in% x & 
                              my.plot.df$degree == "perc_100" & 
                              my.plot.df$variable == "mark" & 
                              my.plot.df$my.mark == type, "value"]))
    
    

  }
  
  
  x <- as.data.frame(x)
  write.table(x, paste0("~/public_html/Borsari_et_al_transdifferentiation_chromatin/ant.del.analysis/perc_25/", type, ".120h.txt"),
              row.names = F, col.names = F, quote=F, sep="\t")
  
}


#********
# BEGIN *
#********

marks <- c("H3K4me3", "H3K4me1", "H3K4me2", "H3K9ac", 
           "H3K27ac", "H3K36me3", "H4K20me1")

my.plot.df <- data.frame(stringsAsFactors = F)

for ( i in marks ){
  
  my.plot.df <- rbind(my.plot.df, my.function2(mark=i))
  
}


my.plot.df$degree <- factor(my.plot.df$degree, levels = c("perc_25", "perc_50", "perc_75", "perc_100"))
my.plot.df$my.mark <- factor(my.plot.df$my.mark, levels = c("H3K27ac", "H3K4me1", "H3K4me2", "H3K9ac",
                                                            "H3K4me3", "H3K36me3", "H4K20me1"))


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.4r.pdf", height=7, width=5.5)
ggplot(my.plot.df[my.plot.df$degree=="perc_25", ], 
       aes(y=my.mark, x=value)) +
  geom_density_ridges(aes(fill = variable), alpha=.5) +
  theme_bw() +
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=13, angle=30, vjust=.5),
        axis.text.y = element_text(size=13),
        strip.text.x = element_text(size=15),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_blank(),
        legend.title = element_blank(),
        legend.text=element_text(size=15),
        legend.position = "top") +
  xlab("time (hours)") +
  geom_vline(xintercept = 12, color="black", linetype ="dashed") +
  scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 72, 120, 168))
dev.off()


# check genes that have 25% upregulation at 120 h
# my.function3(type="expression")
# my.function3(type="H3K4me3")
# my.function3(type="H3K4me1")
# my.function3(type="H3K4me2")
# my.function3(type="H3K27ac")
# my.function3(type="H3K9ac")
# my.function3(type="H3K36me3")
# my.function3(type="H4K20me1")
