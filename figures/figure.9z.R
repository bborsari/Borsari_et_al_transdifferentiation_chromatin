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
  
  df <- read.table(paste0(mark, "/QN.merged/ant.del.analysis/increases/", 
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
marks <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3", "H3K36me3", "H4K20me1")


# 3. list of correspondence between
# degree of upregulation and columns selected
l <- list("25%" = c(1, 5, 9),
          "50%" = c(2, 6, 10),
          "75%" = c(3, 7, 11),
          "100%" = c(4, 8, 12))


# 4. correspondence between 1-12 time-points and 0-168 hours
y <- c("1" = 0, "2" = 3, "3" = 6, "4" = 9,
       "5" = 12, "6" = 18, "7" = 24,
       "8" = 36, "9" = 48, "10" = 72,
       "11" = 120, "12" = 168)


# 5. obtain df for plot - 25% & 100% upregulation
# either considering time-points 1-12
# or hours 0-168
x.25.hours <- data.frame(stringsAsFactors = F)
x.25.tp <- data.frame(stringsAsFactors = F)
x.100.hours <- data.frame(stringsAsFactors = F)
x.100.tp <- data.frame(stringsAsFactors = F)
for ( i in 1:7 ){
  
  x.25.hours <- rbind(x.25.hours, function1(mark=marks[i], degree = "25%", hours = T))
  x.25.tp <- rbind(x.25.tp, function1(mark=marks[i], degree = "25%", hours = F))
  
  x.100.hours <- rbind(x.100.hours, function1(mark=marks[i], degree = "100%", hours = T))
  x.100.tp <- rbind(x.100.tp, function1(mark=marks[i], degree = "100%", hours = F))
  
}


# 6. reorder marks
x.25.hours$mark <- factor(x.25.hours$mark, levels = marks)
x.25.tp$mark <- factor(x.25.tp$mark, levels = marks)
x.100.hours$mark <- factor(x.100.hours$mark, levels = marks)
x.100.tp$mark <- factor(x.100.tp$mark, levels = marks)



# 7. define palette
palette <- c("anticipated" = "#bdbdbd",
             "delayed" = "#403734")



# 8. make plots for 25% up-regulation
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.9z.25.pdf", 
    height=3.5, width=15)
ggplot(x.25.tp, aes(x=group, y=diff, fill=group)) +
  ylab("number of time-points") +
  stat_compare_means(method = "wilcox.test",
                     size = 4.5,
                     label = "p.format",
                     label.x = 1.2)  +  
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
  
