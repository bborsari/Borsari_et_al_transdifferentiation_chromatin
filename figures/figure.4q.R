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


#**********
# PALETTE *
#**********

palette <- c("mark_ant" = "#a1d99b",
             "no_diff" = "#41ab5d",
             "exp_ant" = "#006d2c")




#************
# FUNCTIONS *
#************

hours <- c("1" = 0, "2" = 3, "3" = 6, "4" = 9,
           "5" = 12, "6" = 18, "7" = 24,
           "8" = 36, "9" = 48, "10" = 72,
           "11" = 120, "12" = 168)

my.function1 <- function(degree, expression.mark.m, integer.vector) {
  
  my.list <- list(c(1,2,6), c(1,3,7), c(1,4,8), c(1,5,9))
  names(my.list) <- c("perc_25", "perc_50", "perc_75", "perc_100")
  
  expression.mark.m.subset <- expression.mark.m[, my.list[[degree]]]
  colnames(expression.mark.m.subset) <- c("gene_id", "expression", "mark")
  
  integer.vector.subset <- integer.vector[rownames(integer.vector) %in% expression.mark.m.subset$gene_id, ]
  
  
  exp.ant.df <- expression.mark.m.subset[expression.mark.m.subset$expression < expression.mark.m.subset$mark, ]
  exp.ant.shifts <- c()
  for ( i in 1:nrow(exp.ant.df) ){
    
    my.gene <- exp.ant.df[i, "gene_id"]
    my.exp.tp <- exp.ant.df[i, 2]
    my.alignment <- integer.vector.subset[my.gene, my.exp.tp]
    my.shift <- as.integer(hours[my.exp.tp + my.alignment]) - as.integer(hours[my.exp.tp])
    # my.shift <- as.integer(my.exp.tp + my.alignment) - as.integer(my.exp.tp)
    exp.ant.shifts <- c(exp.ant.shifts, my.shift)
  }
  
  mark.ant.df <- expression.mark.m.subset[expression.mark.m.subset$expression > expression.mark.m.subset$mark, ]
  mark.ant.shifts <- c()
  for ( i in 1:nrow(mark.ant.df) ){
    
    my.gene <- mark.ant.df[i, "gene_id"]
    my.exp.tp <- mark.ant.df[i, 2]
    my.alignment <- integer.vector.subset[my.gene, my.exp.tp]
    my.shift <- as.integer(hours[my.exp.tp + my.alignment]) - as.integer(hours[my.exp.tp])
    # my.shift <- as.integer(my.exp.tp + my.alignment) - as.integer(my.exp.tp)
    mark.ant.shifts <- c(mark.ant.shifts, my.shift)
  }
  
  
  no.diff.df <- expression.mark.m.subset[expression.mark.m.subset$expression == expression.mark.m.subset$mark, ]
  no.diff.shifts <- c()
  for ( i in 1:nrow(no.diff.df) ){
    
    my.gene <- no.diff.df[i, "gene_id"]
    my.exp.tp <- no.diff.df[i, 2]
    my.alignment <- integer.vector.subset[my.gene, my.exp.tp]
    my.shift <- as.integer(hours[my.exp.tp + my.alignment]) - as.integer(hours[my.exp.tp])
    # my.shift <- as.integer(my.exp.tp + my.alignment) - as.integer(my.exp.tp)
    no.diff.shifts <- c(no.diff.shifts, my.shift)
  }
  
  my.out <- data.frame(shift = c(exp.ant.shifts, mark.ant.shifts, no.diff.shifts),
                       group = c(rep("exp_ant", length(exp.ant.shifts)),
                                 rep("mark_ant", length(mark.ant.shifts)),
                                 rep("no_diff", length(no.diff.shifts))),
                       degree = rep(degree, nrow(integer.vector.subset)),
                       gene_id = c(exp.ant.df$gene_id, mark.ant.df$gene_id, no.diff.df$gene_id))
  
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
  
  # 5. read integer vector matrix
  integer.vector <- read.table(paste0(mark, "/QN.merged/ant.del.analysis/", mark, ".integer.vector.tsv"), 
                               h=T, sep = "\t")
  
  # 6. prepare df for plot
  plot.df <- rbind(my.function1(degree="perc_25", 
                                expression.mark.m = expression.mark.m,
                                integer.vector = integer.vector),
                   my.function1(degree="perc_50", 
                                expression.mark.m = expression.mark.m,
                                integer.vector = integer.vector),
                   my.function1(degree="perc_75", 
                                expression.mark.m = expression.mark.m,
                                integer.vector = integer.vector),
                   my.function1(degree="perc_100", 
                                expression.mark.m = expression.mark.m,
                                integer.vector = integer.vector))
  
  plot.df$mark <- mark
  
  return(plot.df)
  
  
  
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


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.4q.pdf", width=4, height = 4)
ggplot(my.plot.df, aes(x=degree, y=shift, color=group)) +
  geom_boxplot() +
  scale_color_manual(values=palette) + 
  theme_bw() +
  theme(axis.title = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_blank()) +
  ylab("NW alignment time lag") +
  xlab("degree of up-regulation") +
  scale_x_discrete(labels = c("25%", "50%", "75%", "100%")) +
  guides(color=F)
dev.off()



