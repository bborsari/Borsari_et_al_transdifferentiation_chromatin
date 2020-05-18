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
  
  df <- read.table(paste0(mark, "/QN.merged/ant.del.analysis/decreases/", 
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
marks <- c("H3K27ac", "H3K9ac", "H4K20me1", "H3K4me3", "H3K4me2", "H3K4me1", "H3K36me3")


# 3. list of correspondence between
# degree of downregulation and columns selected
l <- list("75%" = c(1, 5, 9),
          "50%" = c(2, 6, 10),
          "25%" = c(3, 7, 11),
          "0%" = c(4, 8, 12))


# 4. correspondence between 1-12 time-points and 0-168 hours
hours <- c("1" = 0, "2" = 3, "3" = 6, "4" = 9,
           "5" = 12, "6" = 18, "7" = 24,
           "8" = 36, "9" = 48, "10" = 72,
           "11" = 120, "12" = 168)


# 5. obtain df for plot - 25% & 100% downregulation
x.75 <- data.frame(stringsAsFactors = F)
x.0 <- data.frame(stringsAsFactors = F)
for ( i in 1:7 ){
  
  x.75 <- rbind(x.75, function1(mark=marks[i], degree = "75%"))
  x.0 <- rbind(x.0, function1(mark=marks[i], degree = "0%"))
  
}


# 6. convert time-points to factors
x.75$value <- factor(x.75$value)
x.0$value <- factor(x.0$value)


# 7. define palette
palette <- c("anticipated" = "#bdbdbd",
             "concomitant" = "#737373",
             "delayed" = "#403734")



# 8. make plots for 75% down-regulation

## 8.1. first changing variables (i.e. mark for anticipated, expression for delayed)
x.75.1 <- x.75[x.75$group_variable %in% c("anticipated_mark", "delayed_expression"), ]

### 8.1.1. add p-value of difference to each mark, and reorder marks
sig.75.1 <- c("H3K27ac (ns)", "H3K9ac (ns)", "H4K20me1 (***)",
              "H3K4me3 (ns)", "H3K4me2 (ns)", "H3K4me1 (ns)", "H3K36me3 (ns)")
names(sig.75.1) <- marks
x.75.1$mark <- sig.75.1[x.75.1$mark]
x.75.1$mark <- factor(x.75.1$mark, levels = sig.75.1)
x.75.1$variable <- factor(x.75.1$variable, levels = c("mark", "expression"))


## 8.2. second changing variables (i.e. expression for anticipated, mark for delayed)
x.75.2 <- x.75[x.75$group_variable %in% c("anticipated_expression", "delayed_mark"), ]

### 8.2.1. add p-value of difference to each mark, and reorder marks
sig.75.2 <- c("H3K27ac (ns)", "H3K9ac (ns)", "H4K20me1 (ns)",
              "H3K4me3 (**)", "H3K4me2 (*)", "H3K4me1 (ns)",
              "H3K36me3 (***)")
names(sig.75.2) <- marks
x.75.2$mark <- sig.75.2[x.75.2$mark]
x.75.2$mark <- factor(x.75.2$mark, levels = sig.75.2)
x.75.2$variable <- factor(x.75.2$variable, levels = c("expression", "mark"))



lop <- list()
lop[[1]] <- ggplot(x.75.1, aes(x=value, fill = group)) 
lop[[2]] <- ggplot(x.75.2, aes(x=value, fill = group)) 



lop <- lapply(lop, function(x){x <- x +
  geom_bar() +
  scale_fill_manual(values = palette) +
  theme_bw() +
  theme(axis.title.x = element_text(size=8, hjust = .5),
        axis.title.y = element_text(size=8, vjust = .5),
        axis.text.x = element_text(size=10, angle = 90, vjust=.5),
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
  xlab("time (hours)") +
  ylab("number of genes") +
  scale_x_discrete(drop=FALSE) +
  facet_grid(variable ~ mark) +
  guides(fill=F)})


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.9h.pdf", 
    height=8, width=12)
plot_grid(plotlist = lop[1:2], nrow=2, ncol=1)
dev.off()