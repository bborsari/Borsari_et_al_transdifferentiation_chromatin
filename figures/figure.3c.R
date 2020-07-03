.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/nfs/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


palette <- c("down-\nregulated" = "#810f7c", 
             "bending" = "#737373",
             "up-\nregulated" = "#f16913",
             "peaking" = "#c7e9b4")

# palette2 <- c("positively_correlated" = "#DEA450",
#               "no_peak" = "#9B461F",
#               "negatively_correlated" = "#C9BBC8",
#               "stable" = "#354C3C",
#               "not_correlated" = "#42354C")

palette2 <- c("positively_correlated" = "#5B9F80",
              "no_peak" = "#403734",
              "negatively_correlated" = "#993404",
              "stable" = "#fec44f",
              "not_correlated" = "#ec7014")



#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(reshape2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(rlist)
library(scales)


#********
# BEGIN *
#********

# 1. read marks 6 groups dataframes
marks <- c("H3K27ac", "H3K9ac", "H4K20me1", "H3K4me3", "H3K4me1",
           "H3K36me3", "H3K4me2", "H3K9me3", "H3K27me3")

all.marks.groups <- list()

for ( i in 1:9 ) {
  
  tmp <- read.table(paste0(marks[i], "/QN.merged/", marks[i], ".6.groups.tsv"), h=T, sep="\t", stringsAsFactors = F)
  tmp$group <- gsub("peak_not_TSS", "no_peak", tmp$group)
  tmp$final_class <- gsub("regulation", "-\nregulated", tmp$final_class)
  all.marks.groups[[i]] <- tmp
  
}

rm(tmp)


# # 2. build table of the 6 groups (total set of genes)
# all.marks.groups.table <- data.frame(mark = character(9),
#                                      negatively_correlated = integer(9),
#                                      no_peak = integer(9),
#                                      not_correlated = integer(9),
#                                      positively_correlated = integer(9),
#                                      stable = integer(9),
#                                      stringsAsFactors = F) 
# for (i in 1:9) {
#   tmp <- table(all.marks.groups[[i]]$group)
#   tmp <- c( marks[i], tmp )
#   names(tmp)[1] <- "mark"
#   all.marks.groups.table[i, ] <- tmp
# }



# 3. build table of the 6 groups (divided by upreg/downreg/peaking/bending)
categories <- names(palette)

my.list <- list()

for (cat in categories) {
  
  tmp.df <- data.frame(mark = character(9),
                       negatively_correlated = integer(9),
                       no_peak = integer(9),
                       not_correlated = integer(9),
                       positively_correlated = integer(9),
                       stable = integer(9),
                       stringsAsFactors = F) 
  for (i in 1:9) {
    
    tmp.df2 <- all.marks.groups[[i]]
    tmp.df2$group <- as.factor(tmp.df2$group)
    tmp <- table(tmp.df2[tmp.df2$final_class == cat, "group"])
    tmp <- c( marks[i], tmp )
    names(tmp)[1] <- "mark"
    tmp.df[i, ] <- tmp
    
  } 
  
  tmp.df$class <- cat
  my.list <- list.append(my.list, tmp.df)
  
}


# 4. merge the two tables (the one built on the total set of genes and the one divided by categories)
# all.marks.groups.table$class <- "total"
all.marks.groups.table.class <- rbind(my.list[[1]],
                                      my.list[[2]],
                                      my.list[[3]],
                                      my.list[[4]])

# 5. prepare for plot
all.marks.groups.table.class$negatively_correlated <- as.numeric(all.marks.groups.table.class$negatively_correlated)
all.marks.groups.table.class$no_peak <- as.numeric(all.marks.groups.table.class$no_peak)
all.marks.groups.table.class$not_correlated <- as.numeric(all.marks.groups.table.class$not_correlated)
all.marks.groups.table.class$positively_correlated <- as.numeric(all.marks.groups.table.class$positively_correlated)
all.marks.groups.table.class$stable <- as.numeric(all.marks.groups.table.class$stable)

all.marks.groups.table.class[, 2:6] <- t(apply(all.marks.groups.table.class[, 2:6], 1,
                                               function(x){x/sum(x)}))

all.marks.groups.table.class.melt <- melt(all.marks.groups.table.class)

all.marks.groups.table.class.melt$mark <- factor(all.marks.groups.table.class.melt$mark,
                                                 levels = c("H3K27ac",
                                                            "H3K9ac",
                                                            "H4K20me1",
                                                            "H3K4me3",
                                                            "H3K4me1",
                                                            "H3K36me3",
                                                            "H3K4me2",
                                                            "H3K9me3",
                                                            "H3K27me3"))

all.marks.groups.table.class.melt$class <- factor(all.marks.groups.table.class.melt$class,
                                                  levels = c("up-\nregulated",
                                                             "peaking",
                                                             "down-\nregulated",
                                                             "bending"))

all.marks.groups.table.class.melt$variable <- factor(all.marks.groups.table.class.melt$variable,
                                                     levels = rev(c("no_peak",
                                                                    "positively_correlated",
                                                                    "stable",
                                                                    "not_correlated",
                                                                    "negatively_correlated")))


lop <- list()

for ( i in 1:9 ) {
  lop[[i]] <- ggplot(all.marks.groups.table.class.melt[all.marks.groups.table.class.melt$mark %in% c(marks[i]), ],
                     aes(x=class, y=value, fill=variable)) +
    geom_bar(stat="identity", width=0.8, position = "fill", alpha=1, color="white") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=15),
          axis.title = element_blank(),
          axis.ticks.y = element_line(colour="white"),
          axis.text.y = element_text(size=15, color="white"),
          axis.text.x = element_text(angle=45, hjust=1, size=15),
          strip.text.x = element_text(size=15),
          strip.background.x = element_blank(),
          panel.border = element_rect(color="black"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_fill_manual(values=palette2) + 
    facet_grid(~mark) +
    scale_y_continuous(labels = percent_format()) +
    coord_flip()}


my.legend <- get_legend(lop[[1]])
lop <- lapply(lop, function(x){x <- x+ guides(fill=F)})
lop[[10]] <- my.legend
lop[[1]] <- lop[[1]] + theme(axis.ticks.y = element_line(colour = "black"),
                             axis.text.y = element_text(size=15, colour = "black"))
lop[[6]] <- lop[[6]] + theme(axis.ticks.y = element_line(colour = "black"),
                             axis.text.y = element_text(size=15, colour = "black"))



pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3c.pdf", width=15, height=6.5)
plot_grid(plotlist = lop, nrow=2)
dev.off()
