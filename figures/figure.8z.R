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




function1 <- function(mark, degree) {
  
  df <- read.table(paste0(mark, "/QN.merged/ant.del.analysis/increases/", 
                          "expression.", mark, ".dynamics.tsv"),
                   h=T, sep="\t")
  
  df <- df[, l[[degree]]]
  colnames(df) <- c("expression", "mark", "group")
  
  df.melt <- melt(df, id.vars = "group")
  df.melt$group <- factor(df.melt$group, levels = c("anticipated", "concomitant", "delayed"))
  df.melt$value <- hours[df.melt$value]
  df.melt$group_variable <- ifelse(df.melt$group == "concomitant", 
                                   "concomitant",
                                   paste(df.melt$group, df.melt$variable, sep="_"))

  df.melt$group_variable <- factor(df.melt$group_variable,
                                   levels = c("anticipated_expression",
                                              "delayed_mark",
                                              "concomitant",
                                              "delayed_expression",
                                              "anticipated_mark"))
  df.melt$mark <- mark
  
  return(df.melt)
  
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
hours <- c("1" = 0, "2" = 3, "3" = 6, "4" = 9,
           "5" = 12, "6" = 18, "7" = 24,
           "8" = 36, "9" = 48, "10" = 72,
           "11" = 120, "12" = 168)


# 5. obtain df for plot - 25% & 100% upregulation
x.25 <- data.frame(stringsAsFactors = F)
x.100 <- data.frame(stringsAsFactors = F)
for ( i in 1:7 ){
  
  x.25 <- rbind(x.25, function1(mark=marks[i], degree = "25%"))
  x.100 <- rbind(x.100, function1(mark=marks[i], degree = "100%"))
  
}


# 6. reorder marks
x.25$mark <- factor(x.25$mark, levels = marks)
x.100$mark <- factor(x.100$mark, levels = marks)


# 7. make plots
lop <- list()
lop[[1]] <- ggplot(x.25, aes(x=group_variable, y=value, fill = group_variable)) +
  labs(title = "25% up-regulation")
lop[[2]] <- ggplot(x.100, aes(x=group_variable, y=value, fill = group_variable)) +
  labs(title = "100% up-regulation")

lop <- lapply(lop, function(x){x <- x + 
  geom_violin(scale = "width") +
  theme_bw() +
  coord_flip() +
  theme(axis.title.x = element_text(size=10, hjust = .5),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=10, angle = 30, vjust=.5),
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
        strip.background.x = element_blank(),
        plot.title = element_text(size=13, hjust = .5)) +
  ylab("time (hours)") +
  expand_limits(y=168) +
  scale_y_continuous(labels = c(0, 24, 48, 72, 120, 168),
                     breaks = c(0, 24, 48, 72, 120, 168)) +
  stat_compare_means(comparisons = list(c("anticipated_mark", "delayed_expression"),
                                        c("anticipated_expression", "delayed_mark") ),
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "ns"))) +
  scale_x_discrete(labels = c("expression after mark",
                              "mark after expression",
                              "concomitant",
                              "expression before mark",
                              "mark before expression")) +
  scale_fill_manual(values = c("anticipated_mark" = "blue",
                               "delayed_expression" = "red",
                               "delayed_mark" = "blue",
                               "anticipated_expression" = "red",
                               "concomitant" = "#737373")) +
  facet_grid(~mark) +
  guides(fill=F)})


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8z.pdf", 
    height=7, width=13)
plot_grid(plotlist = lop, nrow=2, ncol=1, align="v")
dev.off()  


