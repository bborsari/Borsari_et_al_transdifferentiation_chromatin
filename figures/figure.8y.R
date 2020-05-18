.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


#************
# LIBRARIES *
#************

library(ggplot2)
library(reshape2)



#********
# BEGIN *
#********

# 1. set working directory
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


# 2. the marks we're analyzing
marks <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3", "H3K36me3", "H4K20me1")


# 3. read df for each mark
x <- data.frame(stringsAsFactors = F)
for ( i in 1:7 ) {
  
  # 3.1. read mark's df
  tmp <- read.table(paste0(marks[i], "/QN.merged/", marks[i], ".6.groups.tsv"),
                    h=T, sep="\t", stringsAsFactors = F)
  
  # 3.2. merge groups "peak_not_TSS" and "no_peak" into "no_peak"   
  tmp$group <- gsub("peak_not_TSS", "no_peak", tmp$group)
  
  # 3.3. keep only upregulated genes
  tmp <- tmp[tmp$final_class == "upregulation", c("group"), drop=F]
  tmp$gene_id <- rownames(tmp)
  
  # 3.4. read mark's dynamics df
  j <- read.table(paste0(marks[i], "/QN.merged/ant.del.analysis2/increases/expression.", 
                         marks[i], ".dynamics.tsv"),
                  h=T, sep="\t", stringsAsFactors = F)
  
  j$gene_id <- rownames(j)
  j <- j[, c("gene_id", "group_25")]
  
  # 3.5. merge tmp and j
  tmp <- merge(tmp, j, by="gene_id", all.x = T)
  tmp$final <- ifelse(is.na(tmp$group_25), tmp$group, tmp$group_25)
  
  
  # 3.6. add mark
  tmp$mark <- marks[i]
  
  
  x <- rbind(x, tmp)

  
}


# 4. reorder groups
x$final <- factor(x$final, levels = c("negatively_correlated",
                                      "not_correlated",
                                      "stable",
                                      "anticipated",
                                      "concomitant",
                                      "delayed",
                                      "no_peak"))


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


# 5.6. subset for upregulated genes
clusters <- clusters[clusters$gene_id %in% tmp$gene_id, ]


# 6. obtain dfs for barplots
all.y <- data.frame(stringsAsFactors = F)

for (k in c("cluster 1", "cluster 2 - 10% least exp.", 
            "cluster 2 - 10% most exp.", "cluster 3")) {
  
  # 6.1. subset m for genes in a specific cluster
  m <- x[x$gene_id %in% clusters[clusters$cluster == k, "gene_id"], ]
  
  # 6.2. compute frequency of each category across marks
  y <- data.frame()
  for ( i in 1:7 ) {
    
    y <- 
      rbind(y, 
            table(m[m$mark == marks[i], "final"]) / length(m[m$mark == marks[i], "final"]) *100)
    
  }
  
  colnames(y) <- levels(x$final)
  rownames(y) <- marks
  
  y <- t(y)
  y.melt <- melt(y)
  y.melt$cluster <- k
  
  all.y <- rbind(all.y, y.melt)
  
}


# 7. reorder clusters
all.y$cluster <- factor(all.y$cluster, levels = c("cluster 3", "cluster 2 - 10% least exp.",
                                                  "cluster 2 - 10% most exp.", "cluster 1"))



# 8. define color palette
palette <- c("no_peak" = "#000000",
             "anticipated" = "#bdbdbd",
             "concomitant" = "#737373",
             "delayed" = "#403734",
             "stable" = "#fec44f",
             "not_correlated" = "#ec7014",
             "negatively_correlated" = "#993404")

# 9. make plot
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8y.pdf", 
    height = 5.5, width = 15)
ggplot(all.y, 
       aes(x=Var2, y=value, fill=Var1)) +
  geom_bar(stat="identity", color = "white", alpha = .9) +
  facet_grid(~cluster) +
  scale_fill_manual(values = palette) +
  theme_bw() +
  theme(axis.title.y = element_text(size=11),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=15, angle = 30, vjust = .5),
        axis.text.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15, angle = 0),
        strip.background = element_blank(),
        plot.title = element_text(size = 22),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  ylab("% of genes")
dev.off()

