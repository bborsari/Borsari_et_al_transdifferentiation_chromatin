.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")



#************
# LIBRARIES *
#************


library(pheatmap)
library(dplyr)
library(tidyr)
library(scales)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(viridis)
library(cowplot)
library(ggplotify)

#************
# FUNCTIONS *
#************

rescale <- function(x){
  return((x -min(x))/(max(x) - min(x)))
}


#********
# BEGIN *
#********

# 1. set working directory
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/")


# 2. read expression matrix of 12448 genes
expression.matrix <- read.table("expression/QN.merged/selected.genes.rep.2.3.after.QN.merged.tsv",
                                h=T, sep="\t")

# 3. retrieve set of not expressed genes
not.expressed.genes <- read.table("expression/silent.genes.txt", h=F, sep="\t", 
                                  stringsAsFactors = F)
not.expressed.genes <- not.expressed.genes$V1


# 4. retrieve set of DE genes
DE.genes <- read.table("expression/QN.merged/expression.matrix.tsv", h=T, sep="\t")
DE.genes <- rownames(DE.genes)


# 5. retrieve set of stably expressed genes 
stably.expressed.genes <- setdiff(rownames(expression.matrix), 
                                  c(not.expressed.genes, DE.genes))

# 6. create a df with gene label (not expressed, stably expressed, DE)
x <- data.frame(gene_id = c(DE.genes, stably.expressed.genes, not.expressed.genes),
                group = c(rep("DE", length(DE.genes)), 
                          rep("stably\ne.", length(stably.expressed.genes)),
                          rep("not\ne.", length(not.expressed.genes))))

# 6. change to HMM wd
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/HMM/marks/")


# 7. analysis of HMM results for the different states
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/HMM/HMM.marks.pdf",
    width = 15, height = 15)
for (i in 2:20) {
  
  # 7.1. read emission matrix
  m1 <- read.table(paste0("HMM.", i, ".response.matrix.tsv"), 
                   h=T, sep="\t")
  
  # 7.2. discard sd values
  m1 <- m1[, paste0("Re", 1:9, "..Intercept.")]
  
  # 7.3. update colnames
  colnames(m1) <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3",
                    "H3K36me3", "H4K20me1", "H3K9me3", "H3K27me3")
  
  # 7.4. sort states according to mean row value
  m1$mean <- apply(m1, 1, mean)
  m1 <- m1[order(m1$mean), ]
  m1$mean <- NULL
  z <- paste0(letters[1:i], ":", gsub("St", "", rownames(m1)))
  names(z) <- rownames(m1)
  
  # 7.5. rescale expression and marks values
  # to range 0-1
  m1 <- as.data.frame(apply(m1, 2, rescale))
  m1$states <- z
  
  palette <- rev(rainbow(i))
  names(palette) <- z
  
  # 7.6. make plot
  p1 <- as.grob(pheatmap(m1[, 1:9], cluster_cols = F, cluster_rows = F,
                         cellheight = 20, cellwidth = 20,
                         silent = T,
                         border_color = NA,
                         main = "emission matrix",
                         color = c('#fff7f3','#fde0dd','#fcc5c0',
                                   '#fa9fb5','#f768a1','#dd3497',
                                   '#ae017e','#7a0177','#49006a'),
                         labels_row = z,
                         annotation_row = m1[, c("states"), drop=F],
                         annotation_colors = list(states = palette))$gtable)
  
  # 7.7. read transition matrix
  m2 <- read.table(paste0("HMM.", i, ".transition.matrix.tsv"), 
                   h=T, sep="\t")
  
  # 7.8. update row- and colnames
  colnames(m2) <- paste0("St", seq(1:i))
  rownames(m2) <- colnames(m2)
  
  # 7.9. reorder states according to emission matrix
  m2 <- m2[rownames(m1), rownames(m1)]
  colnames(m2) <- z
  rownames(m2) <- z
  
  
  # 7.10. make plot
  p2 <- as.grob(pheatmap(m2, cluster_rows = F, cluster_cols = F,
                         cellheight = 20, cellwidth = 20,
                         display_numbers = T,
                         silent = T,
                         number_color = "white",
                         border_color = NA,
                         main = "transition matrix")$gtable)
  
  # 7.11. read matrix of states' sequence for each gene
  m3 <- read.table(paste0("HMM.", i, ".gene.matrix.tsv"), 
                   h=T, sep="\t")
  
  # 7.12. count, for each gene, how many states are found
  m3$states <- apply(m3, 1, function(x){length(table(x))})
  m3$gene_id <- rownames(m3)
  m3 <- merge(m3, x, by = "gene_id")
  m3$group <- factor(m3$group, levels = c("not\ne.", "stably\ne.", "DE"))
  
  y <- melt(table(m3[, c("states", "group")]))
  
  p3 <- ggplot(y, aes(x=group, y=value, label = value, fill = as.factor(states))) +
    geom_bar(stat="identity", position = "fill") +
    labs(fill = "# of states\nin sequence") +
    scale_y_continuous(labels = percent_format()) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 14),
          panel.border = element_rect(color="black"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_fill_brewer()
  
  m3.melt <- melt(m3[, c(1:13, 15)])
  
  m3.melt$value2 <- z[paste0("St", m3.melt$value)]
  
  p4 <- ggplot(m3.melt, aes(x=group, fill = as.factor(value2))) +
    geom_bar(position = "fill") +
    labs(fill = "states") +
    scale_y_continuous(labels = percent_format()) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 14),
          panel.border = element_rect(color="black"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = palette)
  
  
  main.p <- plot_grid(plotlist = list(p1, p2, p3, p4),
                      nrow = 2, scale = c(1, 1, 0.8, 0.8))
  print(main.p)
  
}

dev.off()

