.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


#************
# LIBRARIES *
#************

library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(cowplot)
library(reshape2)



#************
# FUNCTIONS *
#************

rescale <- function(x){
  return((x - min(x))/(max(x) - min(x)))
}


#********
# BEGIN *
#********

# 1. set working directory
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/")


# 2. define color palette
palette <- c("H3K9ac" = "#e7298a",
             "H3K27ac" = "#8e0152",
             "H3K4me3" = "#f1b6da",
             "H3K27me3" = "#253494",
             "H3K9me3" = "#41b6c4",
             "H3K36me3" = "#7fbc41",
             "H4K20me1" = "#276419",
             "H3K4me1" = "#ffbf00",
             "H3K4me2" = "#a67c00",
             "expression" = "black")


# 3. list of marks we're analyzing
marks <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3", "H3K36me3", "H4K20me1")


# 4. read metadata
metadata <- read.table("H3K4me3/QN.merged/metadata.tsv", h=T, sep="\t",
                       stringsAsFactors = F)

# 5. retrieve down-regulated genes
x <- rownames(metadata[metadata$final_class == "downregulation", ])


# 6. read expression matrix & keep only down-regulated genes
expression.m <- read.table("H3K4me3/QN.merged/expression.matrix.tsv", h=T, sep="\t", 
                           stringsAsFactors = F)

# 7. subset expression matrix for down-regulated genes 
expression.m <- expression.m[rownames(expression.m) %in% x, ]


# 8. convert expression profiles to ranges 0-100%
expression.m <- t(apply(expression.m, 1, rescale))


# 9. read marks' matrices
df <- data.frame(stringsAsFactors = F)

for (i in 1:7) {
  
  # 9.1. read mark's matrix
  df.mark <- read.table(paste0(marks[i], "/QN.merged/", marks[i], ".matrix.tsv"),
                        h=T, sep="\t", stringsAsFactors = F)
  
  # 9.2. retrieve genes that are down-regulated and positively correlated 
  # with expression for this mark
  genes <- read.table(paste0(marks[i], "/QN.merged/", marks[i], ".6.groups.tsv"),
                      h=T, sep="\t", stringsAsFactors = F)
  genes <- rownames(genes[genes$final_class == "downregulation" &
                            genes$group == "positively_correlated", ])
  
  # 9.3. subset mark's matrix for the selected genes
  df.mark <- df.mark[rownames(df.mark) %in% genes, ]
  
  # 9.4. rescale mark's matrix in range 0-100%
  df.mark <- t(apply(df.mark, 1, rescale))
  
  # 9.5. subset original expression matrix for the selected genes
  df.expression <- expression.m[rownames(expression.m) %in% genes,]
  
  # 9.6. compute median value in the range 0-100% for mark's matrix
  # for each time-point
  tmp.mark <- as.data.frame(apply(df.mark, 2, median))  
  colnames(tmp.mark) <- "median_value"
  tmp.mark$time_point <- rownames(tmp.mark)
  rownames(tmp.mark) <- 1:12
  tmp.mark$type <- marks[i]
  tmp.mark$type2 <- marks[i]
  df <- rbind(df, tmp.mark)
  
  # 9.7.: repeat 9.6. for expression profiles
  tmp.expression <- as.data.frame(apply(df.expression, 2, median))  
  colnames(tmp.expression) <- "median_value"
  tmp.expression$time_point <- rownames(tmp.expression)
  rownames(tmp.expression) <- 1:12
  tmp.expression$type <- marks[i]
  tmp.expression$type2 <- "expression"
  df <- rbind(df, tmp.expression)
  
  
}

df$type2 <- factor(df$type2, levels = c("expression", "H3K4me1", "H3K4me2", "H3K27ac",
                                        "H3K9ac", "H3K4me3", "H3K36me3", "H4K20me1"))



# 10. make plots
lop <- list()

# 10.1. loop across marks
for (i in 1:7) {
  
  lop[[i]] <- ggplot(df[df$type == marks[i], ], 
                     aes(x=time_point, 
                         y=median_value*100, 
                         group=type2, 
                         color=type2)) +
    geom_point() +
    geom_line() +
    labs(title = marks[i]) +
    theme_bw() +
    scale_color_manual(values = palette) +
    theme(axis.title = element_text(size=11),
          axis.text.x = element_text(size=11, angle = 30, vjust = .5),
          axis.text.y = element_text(size=15),
          strip.text.x = element_text(size=15),
          strip.text.y = element_text(size=15, angle = 0),
          strip.background = element_blank(),
          panel.border = element_rect(color="black"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_blank(),
          plot.title = element_text(size=15, hjust = .5)) +
    scale_x_discrete(labels = c("0", "3", "6", "9", "12", "18",
                                "24", "36", "48", "72", "120", "168")) +
    xlab("time (hours)") +
    guides(color=F, fill=F)
  
  if ( (marks[i] == "H3K4me1") | (marks[i] == "H3K4me3") ) {
    
    lop[[i]] <- lop[[i]] + ylab("% of profile")
    
  } else {
    
    lop[[i]] <- lop[[i]] + ylab("")
    
  }
  
}


# 10.2. store legend
p <- ggplot(df, 
            aes(x=time_point, 
                y=median_value*100, 
                group=type2, 
                color=type2)) +
  geom_point() +
  geom_line() +
  labs(title = marks[i]) +
  theme_bw() +
  scale_color_manual(values = palette) +
  theme(axis.title = element_text(size=11),
        axis.text.x = element_text(size=11, angle = 30, vjust = .5),
        axis.text.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15, angle = 0),
        strip.background = element_blank(),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_blank(),
        plot.title = element_text(size=15, hjust = .5),
        legend.title = element_blank(),
        legend.text = element_text(size=15)) +
  scale_x_discrete(labels = c("0", "3", "6", "9", "12", "18",
                              "24", "36", "48", "72", "120", "168")) +
  xlab("time (hours)")

lop[[8]] <- get_legend(p)



# 10.3. save plot
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8t.pdf",
    width=16, height=8)
plot_grid(plotlist = lop, nrow=2)
dev.off()