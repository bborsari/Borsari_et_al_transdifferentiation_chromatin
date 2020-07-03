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

my.function1 <- function(degree, expression.mark.m, mark) {
  
  my.list <- list(c(1,2,6), c(1,3,7), c(1,4,8), c(1,5,9))
  names(my.list) <- c("perc_75", "perc_50", "perc_25", "perc_0")
  
  expression.mark.m.subset <- expression.mark.m[, my.list[[degree]]]
  colnames(expression.mark.m.subset) <- c("gene_id", "expression", "mark")
  
  group <- c()
  for ( i in 1:nrow(expression.mark.m.subset) ) {
    
    if (expression.mark.m.subset[i, "expression"] < expression.mark.m.subset[i, "mark"]) {
      
      group <- c(group, "exp_ant")
      
    } else if (expression.mark.m.subset[i, "expression"] > expression.mark.m.subset[i, "mark"]) {
      
      group <- c(group, "mark_ant")
      
    } else {
      
      group <- c(group, "no_diff")
    }
    
  }
  
  expression.mark.m.subset$group <- group
  
  expression.mark.m.subset.melt <- melt(expression.mark.m.subset)
  
  
  my.out <- expression.mark.m.subset.melt
  my.out$degree <- degree
  my.out$my.mark <- mark
  
  return(my.out)
  
}


my.function2 <- function(mark, cluster) {
  
  
  # 1. read expression matrix
  expression.m <- read.table(paste0(mark, "/QN.merged/ant.del.analysis/decreases/expression.75.50.25.0.tsv"),
                             h=T, sep="\t")
  expression.m$gene_id <- rownames(expression.m)
  
  # 2. read mark matrix
  mark.m <- read.table(paste0(mark, "/QN.merged/ant.del.analysis/decreases/", mark, ".75.50.25.0.tsv"),
                       h=T, sep="\t")
  mark.m$gene_id <- rownames(mark.m)
  
  # 3. merge mark and expression matrices
  expression.mark.m <- merge(expression.m, mark.m, by="gene_id")
  
  # 4. keep only up-regulated genes
  cluster.genes <- read.table(paste0("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.",
                                     cluster, ".txt"))
  
  metadata <- read.table(paste0(mark, "/QN.merged/metadata.tsv"), 
                         h=T, sep="\t", stringsAsFactors = F)
  
  downreg.genes <- rownames(metadata[metadata$final_class == "downregulation", ])
  downreg.genes <- downreg.genes[downreg.genes %in% cluster.genes$V1]
  
  expression.mark.m <- expression.mark.m[expression.mark.m$gene_id %in% downreg.genes, ]
  
  expression.mark.m[, 2:ncol(expression.mark.m)] <- t(apply(expression.mark.m[, 2:ncol(expression.mark.m)], 1, function(x){x <- hours[x]}))
  
  
  # 6. prepare df for plot
  plot.df <- rbind(my.function1(degree="perc_75", 
                                expression.mark.m = expression.mark.m,
                                mark = mark),
                   my.function1(degree="perc_50", 
                                expression.mark.m = expression.mark.m,
                                mark = mark),
                   my.function1(degree="perc_25", 
                                expression.mark.m = expression.mark.m,
                                mark = mark),
                   my.function1(degree="perc_0", 
                                expression.mark.m = expression.mark.m,
                                mark = mark))
  
  return(plot.df)
  
  
  
}



#********
# BEGIN *
#********

marks <- c("H3K4me3", "H3K4me1", "H3K4me2", "H3K9ac", 
           "H3K27ac", "H3K36me3", "H4K20me1")

my.plot.df.cluster1 <- data.frame(stringsAsFactors = F)
my.plot.df.cluster2 <- data.frame(stringsAsFactors = F)
my.plot.df.cluster3 <- data.frame(stringsAsFactors = F)

for ( i in marks ){
  
  my.plot.df.cluster1 <- rbind(my.plot.df.cluster1, my.function2(mark=i, cluster = 1))
  my.plot.df.cluster2 <- rbind(my.plot.df.cluster2, my.function2(mark=i, cluster = 2))
  my.plot.df.cluster3 <- rbind(my.plot.df.cluster3, my.function2(mark=i, cluster = 3))
  
}


my.plot.df.cluster1$cluster <- "cluster1"
my.plot.df.cluster2$cluster <- "cluster2"
my.plot.df.cluster3$cluster <- "cluster3"

my.plot.df <- rbind(my.plot.df.cluster1,
                    my.plot.df.cluster2,
                    my.plot.df.cluster3)

# plot 75%
df.75 <- my.plot.df[my.plot.df$degree == "perc_75" & my.plot.df$variable == "expression", c("gene_id", "value", "cluster")]
df.75 <- unique(df.75)
df.75$degree <- "75%"

# plot 50%
df.50 <- my.plot.df[my.plot.df$degree == "perc_50" & my.plot.df$variable == "expression", c("gene_id", "value", "cluster")]
df.50 <- unique(df.50)
df.50$degree <- "50%"

# plot 25%
df.25 <- my.plot.df[my.plot.df$degree == "perc_25" & my.plot.df$variable == "expression", c("gene_id", "value", "cluster")]
df.25 <- unique(df.25)
df.25$degree <- "25%"

# plot 0%
df.0 <- my.plot.df[my.plot.df$degree == "perc_0" & my.plot.df$variable == "expression", c("gene_id", "value", "cluster")]
df.0 <- unique(df.0)
df.0$degree <- "0%"


df <- rbind(df.75, df.50, df.25, df.0)
df$degree <- factor(df$degree, levels = c("75%",
                                          "50%",
                                          "25%",
                                          "0%"))
df$cluster <- factor(df$cluster, levels = c("cluster1", "cluster2", "cluster3"))



pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8o.pdf",
    width = 8, height = 4)
ggplot(df, aes(y=value, x=cluster, fill=cluster)) + 
  geom_violin(alpha=.4, colour="white", width = 1) +
  geom_boxplot(width=0.25, alpha = .6, outlier.shape = NA) +
  facet_grid(~degree) +
  theme_bw() +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        strip.text.x = element_text(size=15),
        axis.ticks.x = element_line(color="black"),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_blank(),
        legend.title = element_blank(),
        legend.text=element_text(size=15),
        legend.position = "bottom",
        strip.background.x = element_blank(),
        plot.title = element_text(size=15, hjust = .5)) +
  scale_fill_manual(values = c("#F56E47", "#97BF04", "#772B59")) +
  guides(fill=F) +
  scale_x_discrete(labels = c("1", "2", "3")) +
  xlab("clusters") +
  ylab("time (hours)")
dev.off()

