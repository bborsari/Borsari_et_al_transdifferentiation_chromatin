.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(reshape2)
library(ggalluvial)
library(viridis)


#************
# FUNCTIONS *
#************

function1 <- function(mark) {
  
  df <- read.table(paste0(mark, "/QN.merged/ant.del.analysis/increases/", 
                          "expression.", mark, ".dynamics.tsv"),
                   h=T, sep="\t")
  
  df <- df[, c("group_25", "group_50", "group_75", "group_100")]
  
  df$gene_id <- rownames(df)
  
  df.melt <- melt(df, id.vars = "gene_id")
  df.melt$value <- gsub("concomitant", "co-occurrent", df.melt$value)
  df.melt$value2 <- df.melt$value
  df.melt$value <- factor(df.melt$value, levels = c("anticipated", "co-occurrent", "delayed"))
  df.melt$value2 <- factor(df.melt$value2, levels = c("delayed", "co-occurrent", "anticipated"))
  
  print(mark)
  print(table(df.melt[, c("variable", "value")])/ nrow(df))
  
  p <- ggplot(df.melt,
              aes(x = variable, stratum = value, alluvium = gene_id, label = value)) +
    geom_flow(aes(fill = value2), alpha=0.8) +
    geom_stratum(aes(fill = value), alpha=0.8, color="white") +
    theme_bw() + 
    theme(panel.border = element_rect(color="black"),
          plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_blank(),
          axis.title.y = element_text(size=11),
          axis.text.y = element_text(size=11, angle=90, hjust=.5),
          axis.text.x = element_text(size=12),
          axis.title.x = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size=15),
          legend.title = element_blank(),
          plot.title = element_text(size=17, hjust=.5)) +
    ylab("") +
    scale_x_discrete(labels = c("25%", "50%", "75%", "100%")) +
    scale_color_manual(values=palette) +
    scale_fill_manual(values=palette) +
    labs(title = mark)
  
  return(p)
  
  
}

#********
# BEGIN *
#********

# 1. set working directory
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


# 2. define palette
palette <- c("anticipated" = "#bdbdbd",
             "co-occurrent" = "#737373",
             "delayed" = "#403734")

# 3. the marks we're analyzing
marks <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3", "H3K36me3", "H4K20me1")


# 4. make plots
lop <- list()
for ( i in 1:7 ){
  
  lop[[i]] <- function1(mark=marks[i])
  
}


# 5. add legend
my.legend <- get_legend(lop[[1]])
lop <- lapply(lop, function(x){x <- x+guides(fill=F, color=F)})
lop[[8]] <- my.legend
lop[[1]] <- lop[[1]] + ylab("number of genes")


# 6. save plots
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8v.pdf", 
    width=15, height=3)
p1 <- plot_grid(plotlist = lop[1:7], nrow=1, ncol=7)
plot_grid(plotlist = list(p1, lop[[8]]), nrow=2, rel_heights = c(1, 0.1))
dev.off()



