.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/nfs/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/")

#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)


#********
# BEGIN *
#********

# 1. import dataframe with genes with a profile of cEBPa significantly variable (5% FDR maSigPro)

df <- read.table("cEBPa/QN.merged/cEBPa.QN.merged.maSigPro.out.tsv", h=T, sep="\t")
rownames(df) <- gsub("\\..*", "", rownames(df))

# 2. import cluster dataframes
cl1 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.1.txt", 
                  h=F, sep="\t", stringsAsFactors = F)
cl2 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.2.txt", 
                  h=F, sep="\t", stringsAsFactors = F)
cl3 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.3.txt", 
                  h=F, sep="\t", stringsAsFactors = F)
colnames(cl1) <- "gene_id"
colnames(cl2) <- "gene_id"
colnames(cl3) <- "gene_id"

# 3. proportion of genes with variable profile in the different clusters
cl1.cEBPa <- sum(rownames(df) %in% cl1$gene_id)
cl2.cEBPa <- sum(rownames(df) %in% cl2$gene_id)
cl3.cEBPa <- sum(rownames(df) %in% cl3$gene_id)

pvals <- c()
OR <- c()

# 4. Fisher tests

## cl1 vs cl2
m <- matrix(c(cl1.cEBPa, (nrow(cl1)-cl1.cEBPa), cl2.cEBPa, (nrow(cl2)-cl2.cEBPa)),
            byrow = T, nrow = 2)
pvals <- c(pvals, fisher.test(m)$p.val)
OR <- c(OR, fisher.test(m)$estimate)

## cl2 vs cl3/4
m <- matrix(c(cl2.cEBPa, 
              (nrow(cl2)-cl2.cEBPa), 
              cl3.cEBPa, 
              (nrow(cl3) - cl3.cEBPa)),
            byrow = T, nrow = 2)
pvals <- c(pvals, fisher.test(m)$p.val)
OR <- c(OR, fisher.test(m)$estimate)

## cl1 vs cl3/4
m <- matrix(c(cl1.cEBPa, 
              (nrow(cl1)-cl1.cEBPa), 
              cl3.cEBPa, 
              (nrow(cl3) - cl3.cEBPa)),
            byrow = T, nrow = 2)
pvals <- c(pvals, fisher.test(m)$p.val)
OR <- c(OR, fisher.test(m)$estimate)


plot.df <- data.frame(cluster=c("1", "2", "3"),
                      percentage = c((cl1.cEBPa / nrow(cl1)), 
                                     (cl2.cEBPa / nrow(cl2)),
                                     (cl3.cEBPa / (nrow(cl3)))),
                      stringsAsFactors = F)
plot.df$percentage <- round(as.numeric(plot.df$percentage)*100, 2)

anno <- data.frame(x1 = c(1, 2.1, 1),
                   x2 = c(1.9, 3, 3),
                   y1 = c(27.5, 27.5, 31),
                   y2 = c(28.5, 28.5, 32),
                   cluster = plot.df$cluster,
                   lab = c("3e-10", "1.00", "0.14"),
                   xstar = c(1.45, 2.55, 2),
                   ystar = c(30, 30, 33.5))

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3l.pdf",
    width = 3.5, height=3.5)
ggplot(plot.df, aes(x=cluster, y=percentage, fill=cluster)) + 
  geom_bar(stat = "identity", colour="black") +
  scale_fill_manual(values = c("#F56E47", "#97BF04", "#772B59")) +
  ylab( "% of genes with variable\n cEBPa profile") +
  theme_bw() +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        plot.title = element_text(size=15, hjust = .5),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  guides(fill=F) +
  geom_text(data = anno, aes(x = xstar,  y = ystar, label = lab), size=5) +
  geom_segment(data = anno, aes(x = x1, xend = x1, 
                                y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno, aes(x = x2, xend = x2, 
                                y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno, aes(x = x1, xend = x2, 
                                y = y2, yend = y2),
               colour = "black") +
  xlab("clusters")
dev.off()
