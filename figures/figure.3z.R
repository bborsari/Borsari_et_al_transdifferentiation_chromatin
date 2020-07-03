.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")

palette <- c("positively_correlated" = "#DEA450",
             "no_peak" = "#9B461F",
             "negatively_correlated" = "#a095a0",
             "stable" = "#5d6f62",
             "not_correlated" = "#42354C")



#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(pheatmap)
library(MASS)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(ggrepel)
library(plotly)
library(clustrd)
library(ggalluvial)


#********
# BEGIN *
#********

# 1. read marks 6 groups dataframes
marks <- c("H3K27ac", "H3K9ac", "H4K20me1", "H3K4me3", "H3K4me1",
           "H3K36me3", "H3K4me2", "H3K9me3", "H3K27me3")

all.marks.groups.list <- list()

for ( i in 1:9 ) {
  
  tmp <- read.table(paste0(marks[i], "/QN.merged/", marks[i], ".6.groups.tsv"), h=T, sep="\t", stringsAsFactors = F)
  tmp$group <- gsub("peak_not_TSS", "no_peak", tmp$group)
  tmp$final_class <- gsub("regulation", "-\nregulated", tmp$final_class)
  all.marks.groups.list[[i]] <- tmp
  
}

rm(tmp)


# 2. check order of rownames
for (i in 2:9){
  
  stopifnot(identical(rownames(all.marks.groups.list[[1]]), 
                      rownames(all.marks.groups.list[[i]])))
}


# 3. prepare merged data.frame of groups across all marks
all.marks.groups <- data.frame(H3K27ac = all.marks.groups.list[[1]]$group)

for (i in 2:9) {
  
  tmp <- data.frame(all.marks.groups.list[[i]]$group)
  colnames(tmp) <- marks[i]
  all.marks.groups <- cbind(all.marks.groups, tmp)
  
}

rm(tmp)
rownames(all.marks.groups) <- rownames(all.marks.groups.list[[1]])


# 4. make alluvial plots

# 4.1. - cluster 3
cl3 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.3.txt")
cl3.groups <- all.marks.groups[rownames(all.marks.groups) %in% cl3$V1, ]
cl3.groups$gene_id <- rownames(cl3.groups)
cl3.groups[, 1:9] <- apply(cl3.groups[, 1:9], 2, as.character)
cl3.melt <- melt(cl3.groups, id.vars = "gene_id") 
cl3.melt$value <- factor(cl3.melt$value, levels = c("no_peak", "stable",
                                                    "not_correlated",
                                                    "positively_correlated",
                                                    "negatively_correlated"))

cl3.melt$variable <- factor(cl3.melt$variable, levels = c("H3K9ac", "H3K4me3", "H3K27me3", "H3K27ac",
                                                          "H3K9me3", "H4K20me1", "H3K36me3", "H3K4me2", "H3K4me1"))

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3z.cl3.pdf", width=15, height = 4.5)
ggplot(cl3.melt,
       aes(x = variable, stratum = value, alluvium = gene_id, label = value)) +
  scale_fill_manual(values = palette) +
  scale_color_manual(values =palette) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", aes(color=value, fill=value)) +
  geom_stratum(aes(fill = value)) +
  theme_bw() + 
  labs(title = "cluster 3") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=15),
        axis.text = element_text(size=15),
        axis.title.x = element_blank(),
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        plot.title = element_text(size=15, hjust=.5)) +
  ylab("number of genes")
dev.off()



# 4.2. - cluster 2
cl2 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.2.txt")
cl2.groups <- all.marks.groups[rownames(all.marks.groups) %in% cl2$V1, ]
cl2.groups$gene_id <- rownames(cl2.groups)
cl2.groups[, 1:9] <- apply(cl2.groups[, 1:9], 2, as.character)
cl2.melt <- melt(cl2.groups, id.vars = "gene_id") 
cl2.melt$value <- factor(cl2.melt$value, levels = rev(c("no_peak",
                                                    "negatively_correlated",
                                                    "not_correlated",
                                                    "stable",
                                                    "positively_correlated")))

cl2.melt$value2 <- factor(cl2.melt$value, levels = c("no_peak",
                                                    "negatively_correlated",
                                                    "not_correlated",
                                                    "stable",
                                                    "positively_correlated"))


cl2.melt$variable <- factor(cl2.melt$variable, levels = c("H3K9ac", "H3K27ac", "H3K36me3", 
                                                          "H4K20me1", "H3K4me3", "H3K4me1",
                                                          "H3K4me2","H3K9me3", "H3K27me3"))

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3z.cl2.pdf", width=15, height = 4.5)
ggplot(cl2.melt,
       aes(x = variable, stratum = value, alluvium = gene_id, label = value)) +
  scale_fill_manual(values = palette) +
  scale_color_manual(values =palette) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", aes(color=value2), alpha=.1, fill="white") +
  geom_stratum(aes(fill = value)) +
  theme_bw() + 
  guides(color=F) +
  labs(title = "cluster 2") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=15),
        axis.text = element_text(size=15),
        axis.title.x = element_blank(),
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        plot.title = element_text(size=15, hjust=.5)) +
  ylab("number of genes")
dev.off()


# 4.3. - cluster 1
cl1 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.1.txt")
cl1.groups <- all.marks.groups[rownames(all.marks.groups) %in% cl1$V1, ]
cl1.groups$gene_id <- rownames(cl1.groups)
cl1.groups[, 1:9] <- apply(cl1.groups[, 1:9], 2, as.character)
cl1.melt <- melt(cl1.groups, id.vars = "gene_id") 
cl1.melt$value <- factor(cl1.melt$value, levels = rev(c("positively_correlated",
                                                    "no_peak",
                                                    "negatively_correlated",
                                                    "not_correlated",
                                                    "stable")))

cl1.melt$value2 <- factor(cl1.melt$value, levels = c("positively_correlated",
                                                     "no_peak",
                                                     "negatively_correlated",
                                                     "not_correlated",
                                                     "stable"))


cl1.melt$variable <- factor(cl1.melt$variable, levels = c("H3K36me3", "H3K4me2", "H3K4me3", 
                                                          "H3K4me1", "H3K9ac","H3K27ac", 
                                                          "H4K20me1","H3K9me3", "H3K27me3"))

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3z.cl1.pdf", width=15, height = 4.5)
ggplot(cl1.melt,
       aes(x = variable, stratum = value, alluvium = gene_id, label = value)) +
  scale_fill_manual(values = palette) +
  scale_color_manual(values =palette) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", aes(color=value2), alpha=.1, fill="white") +
  geom_stratum(aes(fill = value)) +
  theme_bw() + 
  guides(color=F) +
  labs(title = "cluster 1") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=15),
        axis.text = element_text(size=15),
        axis.title.x = element_blank(),
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        plot.title = element_text(size=15, hjust=.5)) +
  ylab("number of genes")
dev.off()

