.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


#************
# LIBRARIES *
#************

library(ggplot2)
library(reshape2)
library(ggpubr)



#********
# BEGIN *
#********

# 1. set working directory
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


# 2. read expression matrix
m <- read.table("H3K4me3/QN.merged/expression.matrix.tsv", h=T, sep="\t",
                stringsAsFactors = F)
m$gene_id <- rownames(m)


# 3. read metadata
metadata <- read.table("H3K4me3/QN.merged/metadata.tsv", h=T, sep="\t", 
                       stringsAsFactors = F)


# 4. select up-regulated genes
upreg.genes <- rownames(metadata[metadata$final_class=="upregulation", ])


# 5. read list of genes that belong to:

# 5.1. cluster 3
cluster3 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.3.txt",
                       h=F, sep="\t", stringsAsFactors = F)
cluster3$V2 <- "cluster 3"


# 5.2. cluster 2 lowest
cluster2.lowest <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/genes/cluster2.lowest.txt",
                              h=F, sep="\t", stringsAsFactors = F)
cluster2.lowest$V2 <- "cluster 2 - 10% least exp."


# 5.3. cluster 2 highest
cluster2.highest <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/genes/cluster2.highest.txt",
                               h=F, sep="\t", stringsAsFactors = F)
cluster2.highest$V2 <- "cluster 2 - 10% most exp."


# 5.4. cluster 1
cluster1 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.1.txt",
                       h=F, sep="\t", stringsAsFactors = F)
cluster1$V2 <- "cluster 1"


# 5.5. join all clusters in a unique df
clusters <- rbind(cluster3, cluster2.lowest, cluster2.highest, cluster1)
colnames(clusters) <- c("gene_id", "cluster")


# 6. subset clusters df for up-regulated genes
clusters <- clusters[clusters$gene_id %in% upreg.genes, ]


# 7. join expression and clusters dfs
clusters <- merge(clusters, m[, c("gene_id", "H000", "H168")], by = "gene_id")
clusters$cluster <- factor(clusters$cluster, levels = c("cluster 3", 
                                                        "cluster 2 - 10% least exp.",
                                                        "cluster 2 - 10% most exp.",
                                                        "cluster 1"))

# 8. melt clusters df
clusters.melt <- melt(clusters, id.vars = c("gene_id", "cluster"))


# 9 define color palette
palette <- c("#F56E47", "#cbdf81", "#5a7202", "#772B59")
names(palette) <- c("cluster 1",
                    "cluster 2 - 10% most exp.",
                    "cluster 2 - 10% least exp.",
                    "cluster 3")


# 10. make plot
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.10d.pdf",
    height = 3.5, width = 5)
ggplot(clusters.melt, aes(x=cluster, y = value, fill=cluster)) +
  geom_violin(alpha=.4, colour="white", scale = "width") +
  geom_boxplot(width=0.25, alpha = .6, outlier.shape = NA) +
  stat_compare_means(comparisons = list(c("cluster 3", 
                                          "cluster 2 - 10% least exp."),
                                        c("cluster 2 - 10% most exp.",
                                          "cluster 1"))) +
  ylim(0, 20) +
  facet_grid(.~variable) +
  scale_fill_manual(values = palette) +
  theme_bw() +
  theme(axis.title.x = element_text(size=12),
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
  scale_x_discrete(labels = c("1", "2\n(least)", "2\n(most)", "3")) +
  guides(fill=F) +
  ylab("log2(TPM+1)")
dev.off()
  

