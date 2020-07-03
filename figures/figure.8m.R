.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")

palette <- c("positively_correlated" = "#DEA450",
             "no_peak" = "#9B461F",
             "negatively_correlated" = "#a095a0",
             "stable" = "#5d6f62",
             "not_correlated" = "#42354C")

palette2 <- c("down-regulated" = "#2d7f89", 
              "bending" = "#7acbd5",
              "up-regulated" = "#89372d",
              "peaking" = "#d5847a",
              "flat" = "#8f8f8f")



#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(ggpubr)
library(scales)
library(xtable)



#********
# BEGIN *
#********

# 1. read clusters dataframes
cl1 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.1.txt", h=F, sep="\t")
cl2 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.2.txt", h=F, sep="\t")
cl3 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.3.txt", h=F, sep="\t")

colnames(cl1) <- "gene_id"
colnames(cl2) <- "gene_id"
colnames(cl3) <- "gene_id"


# 2. read metadata file
metadata <- read.table("H3K4me3/QN.merged/metadata.tsv", h=T, sep="\t")
metadata$gene_id <- rownames(metadata)
metadata$final_class <- gsub("regulation", "-regulated", metadata$final_class)


# 3. retrieve number of genes belonging to the different classes within each cluster

## cluster #1
cl1 <- merge(cl1, metadata[, c("gene_id", "final_class")], by = "gene_id")
cl1.v <- c("1", table(cl1$final_class))

## cluster #2
cl2 <- merge(cl2, metadata[, c("gene_id", "final_class")], by = "gene_id")
cl2.v <- c("2", table(cl2$final_class))

## cluster #3
cl3 <- merge(cl3, metadata[, c("gene_id", "final_class")], by = "gene_id")
cl3.v <- c("3", table(cl3$final_class))



# 4. unique df (of % of classes per cluster)
all.clusters <- as.data.frame(rbind(cl1.v, cl2.v, cl3.v))
colnames(all.clusters)[1] <- "cluster"
all.clusters <- apply(all.clusters, 2, as.character)


# 5. fisher test
pvals <- c()
combs <- c()

for (x in c("bending", "down-regulated", "peaking", "up-regulated")) {
  
  # cl1 vs cl2
  m <- matrix(as.numeric(c(all.clusters[1, x], 4995 - as.integer(all.clusters[1, x]),
                           all.clusters[2, x], 2933 - as.integer(all.clusters[2, x]))),
              byrow = T, nrow = 2)
  
  pvals <- c(pvals, fisher.test(m)$p.value)
  combs <- c(paste0(x, "_1vs2"))
  
  # cl2 vs cl3
  m <- matrix(as.numeric(c(all.clusters[2, x], 2933 - as.integer(all.clusters[2, x]),
                           all.clusters[3, x], 102 - as.integer(all.clusters[3, x]))),
              byrow = T, nrow = 2)
  
  pvals <- c(pvals, fisher.test(m)$p.value)
  combs <- c(paste0(x, "_2vs3"))
  
  
  # cl1 vs cl3
  m <- matrix(as.numeric(c(all.clusters[1, x], 4995 - as.integer(all.clusters[1, x]),
                           all.clusters[3, x], 102 - as.integer(all.clusters[3, x]))),
              byrow = T, nrow = 2)
  
  pvals <- c(pvals, fisher.test(m)$p.value)
  combs <- c(paste0(x, "_1vs3"))
  
  
}

all.clusters <- as.data.frame(all.clusters)
all.clusters$cluster <- paste0("cluster", all.clusters$cluster)
all.clusters[, 2:5] <- apply(all.clusters[, 2:5], 2, as.numeric)
cl.n <- c(4995, 2933, 102)
all.clusters[, 2:5] <- apply(all.clusters[, 2:5], 2, function(x){x <- x/cl.n})

all.clusters.melt <- melt(all.clusters, id.vars = "cluster")


pvals <- ifelse(pvals < 0.001, sprintf("%.02e", pvals), round(pvals, 3))

# 6. plots
lop <- list()

# 6.1. bending
anno.bending <- data.frame(x1 = c(1, 2.2, 1),
                           x2 = c(1.8, 3, 3),
                           y1 = c(60, 50, 74),
                           y2 = c(61, 51, 75),
                           cluster = c(1, 2, 3),
                           lab = pvals[1:3],
                           xstar = c(1.4, 2.6, 2),
                           ystar = c(65, 55, 79))

lop[[1]] <- ggplot(all.clusters.melt[all.clusters.melt$variable == "bending", ], 
                   aes(x=cluster, y=value*100)) +
  geom_bar(stat="identity", fill="#7acbd5", color="white") +
  scale_x_discrete(labels = c("1", "2", "3")) +
  geom_text(data = anno.bending, aes(x = xstar,  y = ystar, label = lab), size=4.5) +
  geom_segment(data = anno.bending, aes(x = x1, xend = x1, 
                                   y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno.bending, aes(x = x2, xend = x2, 
                                   y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno.bending, aes(x = x1, xend = x2, 
                                   y = y2, yend = y2),
               colour = "black") +
  theme_bw() + 
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        plot.title = element_text(size = 15, hjust = .5),
        axis.text.y = element_text(size =15),
        axis.text.x = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.background.x = element_blank(),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank()) +
  labs(title = "bending") +
  ylab("% of genes") +
  ylim(0, 80)



# 6.2. down-regulated
anno.downreg <- data.frame(x1 = c(1, 2.2, 1),
                           x2 = c(1.8, 3, 3),
                           y1 = c(60, 50, 74),
                           y2 = c(61, 51, 75),
                           cluster = c(1, 2, 3),
                           lab = pvals[4:6],
                           xstar = c(1.4, 2.6, 2),
                           ystar = c(65, 55, 79))

lop[[2]] <- ggplot(all.clusters.melt[all.clusters.melt$variable == "down-regulated", ], 
                   aes(x=cluster, y=value*100)) +
  geom_bar(stat="identity", fill="#2d7f89", color="white") +
  scale_x_discrete(labels = c("1", "2", "3")) +
  geom_text(data = anno.downreg, aes(x = xstar,  y = ystar, label = lab), size=4.5) +
  geom_segment(data = anno.downreg, aes(x = x1, xend = x1, 
                                        y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno.downreg, aes(x = x2, xend = x2, 
                                        y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno.downreg, aes(x = x1, xend = x2, 
                                        y = y2, yend = y2),
               colour = "black") +
  theme_bw() + 
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        plot.title = element_text(size = 15, hjust = .5),
        axis.text.y = element_text(size =15),
        axis.text.x = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.background.x = element_blank(),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank()) +
  labs(title = "down-regulated") +
  ylab("% of genes") +
  ylim(0, 80)





# 6.3. peaking
anno.peaking <- data.frame(x1 = c(1, 2.2, 1),
                           x2 = c(1.8, 3, 3),
                           y1 = c(60, 50, 74),
                           y2 = c(61, 51, 75),
                           cluster = c(1, 2, 3),
                           lab = pvals[7:9],
                           xstar = c(1.4, 2.6, 2),
                           ystar = c(65, 55, 79))

lop[[3]] <- ggplot(all.clusters.melt[all.clusters.melt$variable == "peaking", ], 
                   aes(x=cluster, y=value*100)) +
  geom_bar(stat="identity", fill="#d5847a", color="white") +
  scale_x_discrete(labels = c("1", "2", "3")) +
  geom_text(data = anno.peaking, aes(x = xstar,  y = ystar, label = lab), size=4.5) +
  geom_segment(data = anno.peaking, aes(x = x1, xend = x1, 
                                        y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno.peaking, aes(x = x2, xend = x2, 
                                        y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno.peaking, aes(x = x1, xend = x2, 
                                        y = y2, yend = y2),
               colour = "black") +
  theme_bw() + 
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        plot.title = element_text(size = 15, hjust = .5),
        axis.text.y = element_text(size =15),
        axis.text.x = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.background.x = element_blank(),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank()) +
  labs(title = "peaking") +
  ylab("% of genes") +
  ylim(0, 80)




# 6.4. up-regulated
anno.upreg <- data.frame(x1 = c(1, 2.2, 1),
                         x2 = c(1.8, 3, 3),
                         y1 = c(60, 50, 74),
                         y2 = c(61, 51, 75),
                         cluster = c(1, 2, 3),
                         lab = pvals[10:12],
                         xstar = c(1.4, 2.6, 2),
                         ystar = c(65, 55, 79))

lop[[4]] <- ggplot(all.clusters.melt[all.clusters.melt$variable == "up-regulated", ], 
                   aes(x=cluster, y=value*100)) +
  geom_bar(stat="identity", fill="#89372d", color="white") +
  scale_x_discrete(labels = c("1", "2", "3")) +
  geom_text(data = anno.upreg, aes(x = xstar,  y = ystar, label = lab), size=4.5) +
  geom_segment(data = anno.upreg, aes(x = x1, xend = x1, 
                                        y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno.upreg, aes(x = x2, xend = x2, 
                                        y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno.upreg, aes(x = x1, xend = x2, 
                                        y = y2, yend = y2),
               colour = "black") +
  theme_bw() + 
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        plot.title = element_text(size = 15, hjust = .5),
        axis.text.y = element_text(size =15),
        axis.text.x = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.background.x = element_blank(),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank()) +
  labs(title = "up-regulated") +
  ylab("% of genes") +
  ylim(0, 80)



pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8m.pdf",
    width = 9, height = 3.5)
plot_grid(plotlist = lop, nrow=1, ncol=4)
dev.off()
