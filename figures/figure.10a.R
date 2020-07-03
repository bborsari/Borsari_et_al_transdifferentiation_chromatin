.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(reshape2)
library(ggpubr)


#************
# FUNCTIONS *
#************




function1 <- function(mark, degree, hours = T) {
  
  df <- read.table(paste0(mark, "/QN.merged/ant.del.analysis/decreases/", 
                          "expression.", mark, ".dynamics.tsv"),
                   h=T, sep="\t")
  
  df <- df[, l[[degree]]]
  colnames(df) <- c("expression", "mark", "group")
  
  if (hours) {
    
    df$expression <- y[df$expression]
    df$mark <- y[df$mark]
    
  }
  
  df$diff <- abs(df$expression - df$mark)
  df <- df[df$group != "concomitant", c("group", "diff")]
  
  df$group <- factor(df$group, levels = c("anticipated", "delayed"))
  
  df$mark <- mark
  
  return(df)
  
}





#********
# BEGIN *
#********

# 1. set working directory
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


# 2. the marks we're analyzing
marks <- c("H3K27ac", "H3K9ac", "H4K20me1", "H3K4me1", "H3K4me3", "H3K36me3", "H3K4me2")


# 3. list of correspondence between
# degree of downregulation and columns selected
l <- list("75%" = c(1, 5, 9),
          "50%" = c(2, 6, 10),
          "25%" = c(3, 7, 11),
          "0%" = c(4, 8, 12))


# 4. correspondence between 1-12 time-points and 0-168 hours
y <- c("1" = 0, "2" = 3, "3" = 6, "4" = 9,
       "5" = 12, "6" = 18, "7" = 24,
       "8" = 36, "9" = 48, "10" = 72,
       "11" = 120, "12" = 168)


# 5. obtain df for plot - 75% & 0% downregulation
# either considering time-points 1-12
# or hours 0-168
x.75.hours <- data.frame(stringsAsFactors = F)
x.75.tp <- data.frame(stringsAsFactors = F)
x.0.hours <- data.frame(stringsAsFactors = F)
x.0.tp <- data.frame(stringsAsFactors = F)
for ( i in 1:7 ){
  
  x.75.hours <- rbind(x.75.hours, function1(mark=marks[i], degree = "75%", hours = T))
  x.75.tp <- rbind(x.75.tp, function1(mark=marks[i], degree = "75%", hours = F))
  
  x.0.hours <- rbind(x.0.hours, function1(mark=marks[i], degree = "0%", hours = T))
  x.0.tp <- rbind(x.0.tp, function1(mark=marks[i], degree = "0%", hours = F))
  
}


# 6. reorder marks
x.75.hours$mark <- factor(x.75.hours$mark, levels = marks)
x.75.tp$mark <- factor(x.75.tp$mark, levels = marks)
x.0.hours$mark <- factor(x.0.hours$mark, levels = marks)
x.0.tp$mark <- factor(x.0.tp$mark, levels = marks)



# 7. define palette
palette <- c("anticipated" = "#bdbdbd",
             "delayed" = "#403734")



# 8. make plots for 75% down-regulation
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.10a.75.pdf", 
    height=3.5, width=15)
ggplot(x.75.tp, aes(x=group, y=diff, fill=group)) +
  ylab("number of time-points") +
  stat_compare_means(method = "wilcox.test",
                     size = 4.5,
                     label = "p.format",
                     label.x = 1.2) +  
  geom_violin(alpha=.6, color = "white") +
  geom_boxplot(width=0.25, alpha = .4, outlier.shape = NA) +
  facet_grid(~mark) +
  scale_fill_manual(values = palette) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=13, vjust = .5),
        axis.text.x = element_text(size=10, vjust=.5),
        axis.text.y = element_text(size=13),
        strip.text.x = element_text(size=15),
        axis.ticks.x = element_line(color="black"),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_blank(),
        legend.title = element_blank(),
        legend.text=element_text(size=15),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size=13),
        plot.title = element_blank()) +
  scale_x_discrete(labels = c("a", "d")) +
  scale_y_continuous(breaks = seq(0, 12, by = 1))
dev.off()

