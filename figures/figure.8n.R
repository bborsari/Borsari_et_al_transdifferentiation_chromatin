.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/ENCODE.data")

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

# 1. list of cell lines
# cell.lines <- read.table("cell.lines.txt", h=F, sep="\t", stringsAsFactors = F)
# cell.lines <- cell.lines$V1
cell.lines <- c("MCF-7", "HepG2", "A549",
                "B.cells", "GM12878", "K562",
                "CD14.positive.monocytes",
                "BlaER", "iMac")
titles <- c("MCF-7", "HepG2", "A549",
            "B cells", "GM12878", "K562",
            "CD14+ monocytes",
            "BlaER (0h)", "iMac (168h)")


# 2. read df with gene quantifications
lom <- list()

for (x in 1:7) {
  
  # 2.1. read rep1
  rep1 <- read.table(paste0(cell.lines[x], 
                            "/quantifications/", 
                            cell.lines[x], 
                            ".1.tsv"), h=T, sep="\t")
  
  # 2.2. read rep2
  rep2 <- read.table(paste0(cell.lines[x], 
                            "/quantifications/", 
                            cell.lines[x], 
                            ".2.tsv"), h=T, sep="\t")
  
  # 2.3. merge replicates
  tmp <- merge(rep1[, c("gene_id", "TPM")], rep2[, c("gene_id", "TPM")],
               by = "gene_id")
  tmp$TPM <- apply(tmp[, 2:3], 1, mean)
  tmp$gene_id <- gsub("\\..*", "", tmp$gene_id)
  rownames(tmp) <- tmp$gene_id
  lom[[x]] <- tmp[, c("gene_id", "TPM"), drop=F]
  
}


# 2.4 add 0h and 168h expression
expression.m <- read.table("../all.marks/expression/QN.merged/expression.matrix.tsv",
                           h=T, sep="\t")

lom[[8]] <- expression.m[, "H000", drop=F]
colnames(lom[[8]]) <- "TPM"
lom[[8]]$gene_id <- rownames(lom[[8]])

lom[[9]] <- expression.m[, "H168", drop=F]
colnames(lom[[9]]) <- "TPM"
lom[[9]]$gene_id <- rownames(lom[[9]])


# 3. read clusters

cl1 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.1.txt")
colnames(cl1) <- "gene_id"
cl1$cluster <- "cluster1"

cl2 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.2.txt")
colnames(cl2) <- "gene_id"
cl2$cluster <- "cluster2"

cl3 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.3.txt")
colnames(cl3) <- "gene_id"
cl3$cluster <- "cluster3"

all.cl <- rbind(cl1, cl2, cl3)


# 4. separate clusters dfs for upreg and downreg

metadata <- read.table("../all.marks/expression/QN.merged/metadata.class2.tsv",
                       h=T, sep="\t")
all.cl.upreg <- all.cl[all.cl$gene_id %in% rownames(metadata[metadata$final_class=="upregulation", ]), ]
all.cl.downreg <- all.cl[all.cl$gene_id %in% rownames(metadata[metadata$final_class=="downregulation", ]), ]


# 5. merge clusters dfs with gene expression dfs
lom.upreg <- list()
lom.downreg <- list()

for ( x in 1:9 ){
  
  lom.upreg[[x]] <- merge(lom[[x]], all.cl.upreg, by = "gene_id")
  lom.downreg[[x]] <- merge(lom[[x]], all.cl.downreg, by = "gene_id")
  
}


# 6. plots
lop.upreg <- list()
lop.downreg <- list()


for ( x in 1:7 ) {
  
  lop.upreg[[x]] <- ggplot(lom.upreg[[x]], aes(x=cluster, y=log2(TPM+1), fill=cluster)) +
                             labs(title = titles[x])
  lop.downreg[[x]] <- ggplot(lom.downreg[[x]], aes(x=cluster, y=log2(TPM+1), fill=cluster)) +
                               labs(title = titles[x])
  
}


for ( x in 8:9 ) {
  
  lop.upreg[[x]] <- ggplot(lom.upreg[[x]], aes(x=cluster, y=TPM, fill=cluster)) +
    labs(title = titles[x])
  lop.downreg[[x]] <- ggplot(lom.downreg[[x]], aes(x=cluster, y=TPM, fill=cluster)) +
    labs(title = titles[x])
  
}


lop.upreg <- lapply(lop.upreg, function(x) {x <- x +
  geom_violin(alpha=.4, colour="white", scale = "width") +
  geom_boxplot(width=0.25, alpha = .6, outlier.shape = NA) +
  stat_compare_means(comparisons = list(c("cluster1", "cluster2"),
                                        c("cluster2", "cluster3"),
                                        c("cluster1", "cluster3"))) +
  scale_x_discrete(labels = c("1", "2", "3")) +
  ylim(0, 20) +
  scale_fill_manual(values = c("#F56E47", "#97BF04", "#772B59")) +
  theme_bw() +
  ylab("log2(TPM+1)") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        plot.title = element_text(size = 15, hjust = .5),
        axis.text.y = element_text(size =12),
        axis.text.x = element_text(size=12),
        strip.text.x = element_text(size=12),
        strip.background.x = element_blank(),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank()) +
  guides(fill=F)} )
  

lop.downreg <- lapply(lop.downreg, function(x) {x <- x +
  geom_violin(alpha=.4, colour="white", scale = "width") +
  geom_boxplot(width=0.25, alpha = .6, outlier.shape = NA) +
  stat_compare_means(comparisons = list(c("cluster1", "cluster2"),
                                        c("cluster2", "cluster3"),
                                        c("cluster1", "cluster3"))) +
  scale_x_discrete(labels = c("1", "2", "3")) +
  ylim(0, 20) +
  scale_fill_manual(values = c("#F56E47", "#97BF04", "#772B59")) +
  theme_bw() + 
  ylab("log2(TPM+1)") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        plot.title = element_text(size = 15, hjust = .5),
        axis.text.y = element_text(size =12),
        axis.text.x = element_text(size=12),
        strip.text.x = element_text(size=12),
        strip.background.x = element_blank(),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank()) +
  guides(fill=F)} )





pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8n.upreg.pdf",
    width = 9, height=9)
plot_grid(plotlist = lop.upreg[c(1,5,6,2,4,7,3,8,9)], nrow=3, ncol=3, align="hv")
dev.off()

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8n.downreg.pdf",
    width = 9, height=9)
plot_grid(plotlist = lop.downreg[c(1,5,6,2,4,7,3,8,9)], nrow=3, ncol=3, align="hv")
dev.off()

