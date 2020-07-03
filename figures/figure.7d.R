.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/nfs/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


#************
# LIBRARIES *
#************

library(ggplot2)
library(ggrepel)


#**********
# PALETTE *
#**********

palette <- c("H3K9ac" = "#af4d85",
             "H3K27ac" = "#630039",
             "H3K4me3" = "#d199b9",
             "H3K27me3" = "#1d2976",
             "H3K9me3" = "#a7add4",
             "H3K36me3" = "#7fbc41",
             "H4K20me1" = "#4c7027",
             "H3K4me1" = "#e5ab00",
             "H3K4me2" = "#a67c00",
             "expression" = "#d9d9d9")


#************
# FUNCTIONS *
#************

my.function <- function(mark=NULL, folder=NULL, expression=F) {
  
  if (expression) {
    
    type = "expression"
    rep1 <- read.table("expression/QN.merged/selected.genes.rep.2.after.QN.merged.tsv", h=T)
    rep2 <- read.table("expression/QN.merged/selected.genes.rep.3.after.QN.merged.tsv", h=T)
    
  } else {
    
    type = mark
    rep1 <- read.table(paste0(folder, "/QN.merged/", mark, ".R1.matrix.after.QN.merged.tsv"), 
                       h=T)
    rep2 <- read.table(paste0(folder, "/QN.merged/", mark, ".R2.matrix.after.QN.merged.tsv"), 
                       h=T)
    
  }
  
  cor.vector <- diag(cor(rep1, rep2))
  cor.df <- data.frame(value = cor.vector,
                       time_point = c(0, 3, 6, 9, 12, 18, 24, 36, 48, 72, 120, 168),
                       type = rep(type, 12))
  return(cor.df)
  
}



#********
# BEGIN *
#********

# 1. obtain Pearson's cc value for expression and the 9 marks
df.plot <- data.frame()

## 1.1 obtain for expression
df.plot <- rbind(df.plot,
                 my.function(expression = T))

## 1.2 obtain for the 9 marks
for (mark in c("H3K9ac", "H3K27ac", "H3K4me3",
               "H3K27me3", "H3K9me3", "H3K36me3",
               "H4K20me1", "H3K4me1", "H3K4me2")) {
  
    df.plot <- rbind(df.plot, my.function(mark = mark, folder = mark))
  
}


# 2. plot
df.plot$pos <- 1
df.plot$time_point <- factor(df.plot$time_point, levels = c(0,3,6,9,12,18,24,36,48,72,120,168))
df.plot$type <- factor(df.plot$type, levels = c("expression",
                                                "H3K27ac",
                                                "H3K9ac",
                                                "H3K4me3",
                                                "H3K36me3",
                                                "H4K20me1",
                                                "H3K4me1",
                                                "H3K4me2",
                                                "H3K9me3",
                                                "H3K27me3"))

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.7d.pdf", width = 14, height = 3)
ggplot(df.plot, aes(x=time_point, y=value, fill=type))  +
  geom_bar(stat="identity", width=0.8) +
  ylab("Pearson's R") +
  xlab("time-points") +
  theme_bw() +
  facet_grid(~type) +
  expand_limits(y=1) +
  theme(panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size=15),
        strip.background.x = element_blank()) +
  guides(fill=F) +
  scale_fill_manual(values = palette)
dev.off()

