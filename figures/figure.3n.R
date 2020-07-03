#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(preprocessCore)


#************
# FUNCTIONS *
#************

palette <- c("H000" = "#cc4c02",
             "H168" = "#005a32")

my.function <- function(mark, category) {
  
      
  setwd(paste0("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/",
               mark, "/QN.merged/aggregation.plots"))
  
  my.df <- read.table(paste0(mark, ".all.groups.agg.plot.tsv"), h=F, sep="\t", quote=NULL)
  colnames(my.df) <- c("position", "value", "group", "rep", "time_point", "type", "class")

  # keep only upreg. or only downreg.
  my.df <- my.df[my.df$class == category, ]
    
  # 1. perform QN normalization across time points and replicates
  my.df2 <- data.frame(H000_rep1 = my.df[my.df$rep == 1 & my.df$time_point == "H000", "value"],
                       H000_rep2 = my.df[my.df$rep == 2 & my.df$time_point == "H000", "value"],
                       H003_rep1 = my.df[my.df$rep == 1 & my.df$time_point == "H003", "value"],
                       H003_rep2 = my.df[my.df$rep == 2 & my.df$time_point == "H003", "value"],
                       H006_rep1 = my.df[my.df$rep == 1 & my.df$time_point == "H006", "value"],
                       H006_rep2 = my.df[my.df$rep == 2 & my.df$time_point == "H006", "value"],
                       H009_rep1 = my.df[my.df$rep == 1 & my.df$time_point == "H009", "value"],
                       H009_rep2 = my.df[my.df$rep == 2 & my.df$time_point == "H009", "value"],
                       H012_rep1 = my.df[my.df$rep == 1 & my.df$time_point == "H012", "value"],
                       H012_rep2 = my.df[my.df$rep == 2 & my.df$time_point == "H012", "value"],
                       H018_rep1 = my.df[my.df$rep == 1 & my.df$time_point == "H018", "value"],
                       H018_rep2 = my.df[my.df$rep == 2 & my.df$time_point == "H018", "value"],
                       H024_rep1 = my.df[my.df$rep == 1 & my.df$time_point == "H024", "value"],
                       H024_rep2 = my.df[my.df$rep == 2 & my.df$time_point == "H024", "value"],
                       H036_rep1 = my.df[my.df$rep == 1 & my.df$time_point == "H036", "value"],
                       H036_rep2 = my.df[my.df$rep == 2 & my.df$time_point == "H036", "value"],
                       H048_rep1 = my.df[my.df$rep == 1 & my.df$time_point == "H048", "value"],
                       H048_rep2 = my.df[my.df$rep == 2 & my.df$time_point == "H048", "value"],
                       H072_rep1 = my.df[my.df$rep == 1 & my.df$time_point == "H072", "value"],
                       H072_rep2 = my.df[my.df$rep == 2 & my.df$time_point == "H072", "value"],
                       H120_rep1 = my.df[my.df$rep == 1 & my.df$time_point == "H120", "value"],
                       H120_rep2 = my.df[my.df$rep == 2 & my.df$time_point == "H120", "value"],
                       H168_rep1 = my.df[my.df$rep == 1 & my.df$time_point == "H168", "value"],
                       H168_rep2 = my.df[my.df$rep == 2 & my.df$time_point == "H168", "value"])
  
  my.df2 <- as.matrix(my.df2)
  my.df2.QN  <- normalize.quantiles(my.df2, copy = TRUE)
  rownames(my.df2.QN) <- rownames(my.df2)
  colnames(my.df2.QN) <- colnames(my.df2)
  my.df2.QN <- as.data.frame(my.df2.QN)
  
  # compute mean of H000 rep 1&2, H168 rep 1&2 after QN
  my.df2.QN <- my.df2.QN[, c("H000_rep1", "H000_rep2", "H168_rep1", "H168_rep2")]
  my.df2.QN$H000_mean <- (my.df2.QN$H000_rep1 + my.df2.QN$H000_rep2) /2
  my.df2.QN$H168_mean <- (my.df2.QN$H168_rep1 + my.df2.QN$H168_rep2) /2
  
  # keep only rep.1
  my.df <- my.df[my.df$rep == 1 & my.df$time_point %in% c("H000", "H168"), ]
  
  # substitute rep.1 values with mean across replicates
  my.df[my.df$time_point == "H000", "value"] <- my.df2.QN$H000_mean
  my.df[my.df$time_point == "H168", "value"] <- my.df2.QN$H168_mean
  
  # keep only major TSS
  my.df <- my.df[my.df$type == "major", ]
  
  return(my.df)
  
}


my.function2 <- function(my.df, mark){
  
  my.x.lim = c(-2000, +2000)
  my.x.intercept = c(-250, 500)
  
  groups <- c("positively_correlated", "negatively_correlated", "not_correlated",
              "stable", "no_peak")
  
  this.list <- list()
  
  i=1
  
  for (k in groups) {
    
    if (nrow(my.df[my.df$group == k, ]) > 0) {
      
      
      this.list[[i]] <- ggplot(my.df[my.df$group==k, ],
                               aes(x=position, y=value, group=time_point, colour=time_point)) 
      i <- i+1
      
    }
    
  }

    
  my.y <- max(as.numeric(my.df[my.df$time_point  %in% c("H000", "H168"), "value"]), na.rm = T)
  
  
  this.list <- lapply(this.list, function(x) { x <- x + facet_grid(. ~ group) + geom_line() +
    guides(colour=F) +
    scale_color_manual(values=palette) +
    scale_x_continuous(limits = my.x.lim) +
    scale_y_continuous(limits = c(0, my.y+0.1)) +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(size=10, angle = 30, vjust = 0.5),
          axis.text.y = element_text(size = 10)) +
    geom_vline(xintercept = my.x.intercept, linetype = "dashed")}) 
  
  return(this.list)
  
}


#********************
# AGGREGATION PLOTS *
#********************


###------------------
### retrieve legend
###------------------

legend.df <- data.frame(time_point = c("H000", "H168"))
p <- ggplot(legend.df, aes(x=time_point, colour =time_point)) + 
  geom_density() + scale_color_manual(values=palette, labels = c("0h", "168h")) +
  labs(colour = "Time point")

my.legend <- ggdraw(get_legend(p))

###----------
### H3K27ac
###----------

# 1. upregulated genes
H3K27ac.mark.upreg <- my.function(mark = "H3K27ac", category = "upregulation")
p.H3K27ac.mark.upreg <- plot_grid(plotlist = my.function2(my.df = H3K27ac.mark.upreg,
                                                          mark = "H3K27ac"),
                                  nrow=1, ncol=5)


# 2. downregulated genes
H3K27ac.mark.downreg <- my.function(mark = "H3K27ac", category = "downregulation")
p.H3K27ac.mark.downreg <- plot_grid(plotlist = my.function2(my.df = H3K27ac.mark.downreg,
                                                            mark = "H3K27ac"),
                                    nrow=1, ncol=5)


# 3. upregulated + downregulated genes
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3n.H3K27ac.pdf", 
    width=16, height=8)
plot_grid(plotlist = list(p.H3K27ac.mark.upreg, p.H3K27ac.mark.downreg),
          nrow=2, ncol=1, scale = .9, labels = c("upregulated genes", "downregulated genes"))
dev.off()




###----------
### H3K9ac
###----------

# 1. upregulated genes
H3K9ac.mark.upreg <- my.function(mark = "H3K9ac", category = "upregulation")
p.H3K9ac.mark.upreg <- plot_grid(plotlist = my.function2(my.df = H3K9ac.mark.upreg,
                                                          mark = "H3K9ac"),
                                  nrow=1, ncol=5)


# 2. downregulated genes
H3K9ac.mark.downreg <- my.function(mark = "H3K9ac", category = "downregulation")
p.H3K9ac.mark.downreg <- plot_grid(plotlist = my.function2(my.df = H3K9ac.mark.downreg,
                                                            mark = "H3K9ac"),
                                    nrow=1, ncol=5)


# 3. upregulated + downregulated genes
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3n.H3K9ac.pdf", 
    width=16, height=8)
plot_grid(plotlist = list(p.H3K9ac.mark.upreg, p.H3K9ac.mark.downreg),
          nrow=2, ncol=1, scale = .9, labels = c("upregulated genes", "downregulated genes"))
dev.off()




###----------
### H3K4me3
###----------

# 1. upregulated genes
H3K4me3.mark.upreg <- my.function(mark = "H3K4me3", category = "upregulation")
p.H3K4me3.mark.upreg <- plot_grid(plotlist = my.function2(my.df = H3K4me3.mark.upreg,
                                                          mark = "H3K4me3"),
                                  nrow=1, ncol=5)


# 2. downregulated genes
H3K4me3.mark.downreg <- my.function(mark = "H3K4me3", category = "downregulation")
p.H3K4me3.mark.downreg <- plot_grid(plotlist = my.function2(my.df = H3K4me3.mark.downreg,
                                                            mark = "H3K4me3"),
                                    nrow=1, ncol=5)


# 3. upregulated + downregulated genes
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3n.H3K4me3.pdf", 
    width=16, height=8)
plot_grid(plotlist = list(p.H3K4me3.mark.upreg, p.H3K4me3.mark.downreg),
          nrow=2, ncol=1, scale = .9, labels = c("upregulated genes", "downregulated genes"))
dev.off()




###----------
### H3K4me1
###----------

# 1. upregulated genes
H3K4me1.mark.upreg <- my.function(mark = "H3K4me1", category = "upregulation")
p.H3K4me1.mark.upreg <- plot_grid(plotlist = my.function2(my.df = H3K4me1.mark.upreg,
                                                          mark = "H3K4me1"),
                                  nrow=1, ncol=5)


# 2. downregulated genes
H3K4me1.mark.downreg <- my.function(mark = "H3K4me1", category = "downregulation")
p.H3K4me1.mark.downreg <- plot_grid(plotlist = my.function2(my.df = H3K4me1.mark.downreg,
                                                            mark = "H3K4me1"),
                                    nrow=1, ncol=5)


# 3. upregulated + downregulated genes
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3n.H3K4me1.pdf", 
    width=16, height=8)
plot_grid(plotlist = list(p.H3K4me1.mark.upreg, p.H3K4me1.mark.downreg),
          nrow=2, ncol=1, scale = .9, labels = c("upregulated genes", "downregulated genes"))
dev.off()




###----------
### H3K4me2
###----------

# 1. upregulated genes
H3K4me2.mark.upreg <- my.function(mark = "H3K4me2", category = "upregulation")
p.H3K4me2.mark.upreg <- plot_grid(plotlist = my.function2(my.df = H3K4me2.mark.upreg,
                                                          mark = "H3K4me2"),
                                  nrow=1, ncol=5)


# 2. downregulated genes
H3K4me2.mark.downreg <- my.function(mark = "H3K4me2", category = "downregulation")
p.H3K4me2.mark.downreg <- plot_grid(plotlist = my.function2(my.df = H3K4me2.mark.downreg,
                                                            mark = "H3K4me2"),
                                    nrow=1, ncol=5)


# 3. upregulated + downregulated genes
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3n.H3K4me2.pdf", 
    width=16, height=8)
plot_grid(plotlist = list(p.H3K4me2.mark.upreg, p.H3K4me2.mark.downreg),
          nrow=2, ncol=1, scale = .9, labels = c("upregulated genes", "downregulated genes"))
dev.off()




###----------
### H3K27me3
###----------

# 1. upregulated genes
H3K27me3.mark.upreg <- my.function(mark = "H3K27me3", category = "upregulation")
p.H3K27me3.mark.upreg <- plot_grid(plotlist = my.function2(my.df = H3K27me3.mark.upreg,
                                                          mark = "H3K27me3"),
                                  nrow=1, ncol=5)


# 2. downregulated genes
H3K27me3.mark.downreg <- my.function(mark = "H3K27me3", category = "downregulation")
p.H3K27me3.mark.downreg <- plot_grid(plotlist = my.function2(my.df = H3K27me3.mark.downreg,
                                                            mark = "H3K27me3"),
                                    nrow=1, ncol=5)


# 3. upregulated + downregulated genes
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3n.H3K27me3.pdf", 
    width=16, height=8)
plot_grid(plotlist = list(p.H3K27me3.mark.upreg, p.H3K27me3.mark.downreg),
          nrow=2, ncol=1, scale = .9, labels = c("upregulated genes", "downregulated genes"))
dev.off()





###----------
### H3K9me3
###----------

# 1. upregulated genes
H3K9me3.mark.upreg <- my.function(mark = "H3K9me3", category = "upregulation")
p.H3K9me3.mark.upreg <- plot_grid(plotlist = my.function2(my.df = H3K9me3.mark.upreg,
                                                          mark = "H3K9me3"),
                                  nrow=1, ncol=5)


# 2. downregulated genes
H3K9me3.mark.downreg <- my.function(mark = "H3K9me3", category = "downregulation")
p.H3K9me3.mark.downreg <- plot_grid(plotlist = my.function2(my.df = H3K9me3.mark.downreg,
                                                            mark = "H3K9me3"),
                                    nrow=1, ncol=5)


# 3. upregulated + downregulated genes
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3n.H3K9me3.pdf", 
    width=16, height=8)
plot_grid(plotlist = list(p.H3K9me3.mark.upreg, p.H3K9me3.mark.downreg),
          nrow=2, ncol=1, scale = .9, labels = c("upregulated genes", "downregulated genes"))
dev.off()


