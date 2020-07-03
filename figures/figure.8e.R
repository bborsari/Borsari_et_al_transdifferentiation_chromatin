.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/nfs/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/Hi-C/")

#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(pheatmap)
library(reshape2)
library(dplyr)


#************
# FUNCTIONS *
#************

# this function is to retrieve the % of
# genes with at least 1 A/B transition 
# in the different clusters
my.function1 <- function(f1.df, f1.title) {
  
  # 1. keep only columns corresponding to gene_id, number of transitions and cluster info
  f1.df <- f1.df[, c("gene_id", "number_transitions", "cluster")]
  
  # 2. sort the df so that for genes with multiple entries we 
  # prioritize the entry with highest number of switches
  f1.df <- f1.df[order(f1.df$gene_id,
                       f1.df$number_transitions,
                       decreasing = T), ]
  
  # 3. for genes with multiple entries,
  # keep the one with highest number of switches
  f1.df <- f1.df %>%
    group_by(gene_id) %>%
    top_n(wt = number_transitions, n = 1) %>%
    unique()
  
  # 4. count at most one transition per genes
  f1.df$number_transitions <- ifelse(f1.df$number_transitions > 0, 1, 0)
  
  # 5. summarize the genes with and w/o transitions per cluster
  print(table(f1.df[, c("number_transitions", "cluster")]))
  
  
  f1.pvals <- c()
  
  # 6. count the numbers of genes with (x) 
  # and without (y) transitions per cluster
  
  ## cluster 1
  f1.x1 <- nrow(f1.df[f1.df$cluster=="cluster1" & f1.df$number_transitions>0, ])
  f1.y1 <- sum(f1.df$cluster=="cluster1") - f1.x1
  
  ## cluster 2
  f1.x2 <- nrow(f1.df[f1.df$cluster=="cluster2" & f1.df$number_transitions>0, ])
  f1.y2 <- sum(f1.df$cluster=="cluster2") - f1.x2
  
  ## cluster 3
  f1.x3 <- nrow(f1.df[f1.df$cluster=="cluster3" & f1.df$number_transitions>0, ])
  f1.y3 <- sum(f1.df$cluster=="cluster3") - f1.x3
  
  
  # 7. comparisons of the % of genes
  # with at least one transition between clusters
  
  ## cluster 1 vs. 2
  f1.m.1.2 <- matrix(c(f1.x1, f1.y1, f1.x2, f1.y2), nrow = 2, byrow = T)
  f1.pvals <- c(f1.pvals, fisher.test(f1.m.1.2)$p.value)
  
  f1.plot.df.1.2 <- data.frame(cluster=c("1", "2"),
                               percentage = c(f1.x1/(f1.x1+f1.y1), 
                                              f1.x2/(f1.x2+f1.y2)),
                               stringsAsFactors = F)
  
  ## cluster 2 vs. 3
  f1.m.2.3 <- matrix(c(f1.x2, f1.y2, f1.x3, f1.y3), nrow = 2, byrow = T)
  f1.pvals <- c(f1.pvals, fisher.test(f1.m.2.3)$p.value)
  
  f1.plot.df.2.3 <- data.frame(cluster=c("2", "3"),
                               percentage = c(f1.x2/(f1.x2+f1.y2), 
                                              f1.x3/(f1.x3+f1.y3)),
                               stringsAsFactors = F)
  
  ## cluster 1 vs. 3
  f1.m.1.3 <- matrix(c(f1.x1, f1.y1, f1.x3, f1.y3), nrow = 2, byrow = T)
  f1.pvals <- c(f1.pvals, fisher.test(f1.m.1.3)$p.value)
  
  f1.plot.df.1.3 <- data.frame(cluster=c("1", "3"),
                               percentage = c(f1.x1/(f1.x1+f1.y1), 
                                              f1.x3/(f1.x3+f1.y3)),
                               stringsAsFactors = F)
  
  
  # 8. make plot
  
  f1.plot.df <- rbind(f1.plot.df.1.2,
                      f1.plot.df.2.3,
                      f1.plot.df.1.3)
  
  f1.plot.df <- unique(f1.plot.df)
  
  f1.plot.df$percentage <- round(as.numeric(f1.plot.df$percentage)*100, 2)
  
  f1.anno <- data.frame(x1 = c(1, 2, 1),
                        x2 = c(2, 3, 3),
                        y1 = c(18, 27, 31),
                        y2 = c(19, 28, 32),
                        cluster = f1.plot.df$cluster,
                        lab = round(f1.pvals, 3),
                        xstar = c(1.5, 2.5, 2),
                        ystar = c(21, 30, 34))
  
  print(f1.plot.df)
  
  f1.p <- ggplot(f1.plot.df, aes(x=cluster, y=percentage, fill=cluster)) + 
    geom_bar(stat = "identity", colour="black") +
    ylab( "% of genes with at least \n1 A/B transition") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=15),
          axis.text = element_text(size=15),
          plot.title = element_text(size=15, hjust = .5),
          panel.border = element_rect(color="black"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_fill_manual(values = c("#F56E47", "#97BF04", "#772B59")) +
    labs(title = f1.title) +
    guides(fill=F) +
    geom_text(data = f1.anno, aes(x = xstar,  y = ystar, label = lab), size=5) +
    geom_segment(data = f1.anno, aes(x = x1, xend = x1, 
                                  y = y1, yend = y2),
                 colour = "black") +
    geom_segment(data = f1.anno, aes(x = x2, xend = x2, 
                                  y = y1, yend = y2),
                 colour = "black") +
    geom_segment(data = f1.anno, aes(x = x1, xend = x2, 
                                  y = y2, yend = y2),
                 colour = "black")
  
  return(f1.p)
  
}



#********
# BEGIN *
#********

# 1. import complete datasets with A/B values for all time-points
# + number of transitions
all.tp <- read.table("all.time.points/complete.data.set.tsv", h=F, sep="\t")
matched.tp <- read.table("matching.time.points/complete.data.set.tsv", h=F, sep="\t")


# 2. update colnames
colnames(all.tp) <- c("gene_id", "0h", "9h",
                     "18h", "24h", "48h",
                     "72h", "96h", "120h",
                     "144h", "168h", "192h", "number_transitions")

colnames(matched.tp) <- c("gene_id", "0h", "9h",
                          "18h", "24h", "48h",
                          "72h", "120h",
                          "168h", "number_transitions")



# 3. import list of genes per cluster
clusters <- read.table("cluster.1.2.3.txt")

# 4. update colnames
colnames(clusters) <- c("gene_id", "cluster")

# 5. add to table of A/B values + transitions
# the info regarding the clusters
all.tp <- merge(all.tp, clusters, all.x = T, by = "gene_id")
matched.tp <- merge(matched.tp, clusters, all.x = T, by = "gene_id")


# 6. % of genes with at least 1 transition
lop <- list()
lop[[1]] <- my.function1(f1.df = all.tp, f1.title = "all time-points")
lop[[2]] <- my.function1(f1.df = matched.tp, f1.title = "matched time-points")

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8e.pdf",
    width = 3.5, height=3.5)
#plot_grid(plotlist = lop, nrow=1, ncol=2, align="h")
lop[[1]] <- lop[[1]] + 
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size=15)) + 
  xlab("clusters")
lop[[1]]
dev.off()








