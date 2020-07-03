.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")

palette <- c("positively_correlated" = "#5B9F80",
             "no_peak" = "#403734",
             "negatively_correlated" = "#993404",
             "stable" = "#fec44f",
             "not_correlated" = "#ec7014")



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
library(ggExtra)
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


# 4. cluster correspondence analysis
# performed in the cluster with R-3.5.2
# with the command below


# out.mca <- clusmca(all.marks.groups, nclus = 3, ndim = 3, method=c("clusCA","iFCB","MCAk"),
#         alphak = .5, nstart = 100, smartStart = NULL, gamma = TRUE,
#         seed = 1234)


load(".mca.RData")



# 5. prepare legend
cl3 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.3.txt")
cl3.groups <- all.marks.groups[rownames(all.marks.groups) %in% cl3$V1, ]
cl3.groups$gene_id <- rownames(cl3.groups)
cl3.groups[, 1:9] <- apply(cl3.groups[, 1:9], 2, as.character)
cl3.melt <- melt(cl3.groups, id.vars = "gene_id") 
cl3.melt$value <- factor(cl3.melt$value, levels = c("negatively_correlated",
                                                    "not_correlated",
                                                    "stable",
                                                    "positively_correlated",
                                                    "no_peak"
                                                    ))

cl3.melt$variable <- factor(cl3.melt$variable, levels = c("H3K9ac", "H3K4me3", "H3K27me3", "H3K27ac",
                                                          "H3K9me3", "H4K20me1", "H3K36me3", "H3K4me2", "H3K4me1"))

palette <- c("positively_correlated" = "#5B9F80",
             "no_peak" = "#403734",
             "negatively_correlated" = "#993404",
             "stable" = "#fec44f",
             "not_correlated" = "#ec7014")


plot.legend <- ggplot(cl3.melt,
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
        plot.title = element_text(size=15, hjust=.5),
        legend.position = "bottom") +
  ylab("number of genes")

my.legend <- get_legend(plot.legend)



# 6. plot PC1 vs. PC2
## 6.1. 3d plot of the variables (marks)
cats = apply(all.marks.groups, 2, function(x) nlevels(as.factor(x)))
mca_vars_df = data.frame(out.mca$attcoord, mark = rep(names(cats), cats))
mca_vars_df$group <- c(names(table(all.marks.groups$H3K27ac)),
                       names(table(all.marks.groups$H3K9ac)),
                       names(table(all.marks.groups$H4K20me1)),
                       names(table(all.marks.groups$H3K4me3)),
                       names(table(all.marks.groups$H3K4me1)),
                       names(table(all.marks.groups$H3K36me3)),
                       names(table(all.marks.groups$H3K4me2)),
                       names(table(all.marks.groups$H3K9me3)),
                       names(table(all.marks.groups$H3K27me3)))
rownames(mca_vars_df) <- paste(mca_vars_df$mark, mca_vars_df$group, sep=".")
mca_vars_df$frequency <- c(round(table(all.marks.groups$H3K27ac) / nrow(all.marks.groups), 3)*100,
                           round(table(all.marks.groups$H3K9ac) / nrow(all.marks.groups), 3)*100,
                           round(table(all.marks.groups$H4K20me1) / nrow(all.marks.groups), 3)*100,
                           round(table(all.marks.groups$H3K4me3) / nrow(all.marks.groups), 3)*100,
                           round(table(all.marks.groups$H3K4me1) / nrow(all.marks.groups), 3)*100,
                           round(table(all.marks.groups$H3K36me3) / nrow(all.marks.groups), 3)*100,
                           round(table(all.marks.groups$H3K4me2) / nrow(all.marks.groups), 3)*100,
                           round(table(all.marks.groups$H3K9me3) / nrow(all.marks.groups), 3)*100,
                           round(table(all.marks.groups$H3K27me3) / nrow(all.marks.groups), 3)*100)
colnames(mca_vars_df)[1:3] <- c("X", "Y", "Z")

mca_vars_df$alpha <- mca_vars_df$frequency
mca_vars_df[mca_vars_df$mark == "H3K27me3" & mca_vars_df$group=="negatively_correlated", "alpha"] <- 40 
mca_vars_df[mca_vars_df$mark %in% c("H3K27ac",
                                    "H3K9ac",
                                    "H3K4me1",
                                    "H3K4me2",
                                    "H3K4me3") & mca_vars_df$group=="no_peak", "alpha"] <- 30 


# p1 <- ggplot(mca_vars_df, aes(x=X, y=Y, fill=group, label=mark)) +
#   geom_point(aes(size=frequency), color="black",shape=21) +
#   scale_color_manual(values = palette) +
#   scale_fill_manual(values = palette) +
#   geom_text_repel(aes(alpha=alpha, color=group), size=4.5) +
#   theme_bw() + 
#   coord_cartesian() +
#   guides(alpha=F, color=F, fill=F) +
#   xlim(-5, 15) +
#   ylim(-1.8, 2.5) +
#   theme(panel.border = element_rect(color="black"), 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         axis.title = element_text(size=15),
#         axis.text = element_text(size=15),
#         legend.text = element_text(size=15),
#         legend.title = element_text(size=15),
#         legend.position = c(.9, .8))
palette2 <- c("cluster1" = "#F56E47", "cluster2" = "#97BF04", "cluster3" = "#772B59")

test <- data.frame(X = mca_vars_df$X,
                   Z = -mca_vars_df$Z,
                   Y = -mca_vars_df$Y)
shapes <- c("H3K4me1" = 0,
            "H3K4me2" = 7,
            "H3K27ac" = 1,
            "H3K9ac" = 10,
            "H3K4me3" = 2,
            "H3K36me3" = 5,
            "H4K20me1" = 9,
            "H3K9me3" = 4,
            "H3K27me3" = 8)

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3g.marks.pdf", 
    width=5, height=4)
scatterplot3d(test,
              pch=shapes[mca_vars_df$mark],
              color=palette[mca_vars_df$group],
              angle = 50, grid = F, box = T,
              cex.symbols = 0.8,
              scale.y = 0.5)
dev.off()



## 6.2. 3d plot of the genes
mca_genes_df  <- as.data.frame(out.mca$obscoord)
colnames(mca_genes_df) <- c("X", "Y", "Z")
rownames(mca_genes_df) <- rownames(all.marks.groups)
mca_genes_df$cluster <- out.mca$cluster
mca_genes_df$cluster <- as.factor(mca_genes_df$cluster)
mca_genes_df$cluster <- paste0("cluster", mca_genes_df$cluster)


# p2 <- ggplot(mca_genes_df, aes(x=X, y=Y, color=cluster, fill=cluster)) +
#   stat_density_2d() +
#   scale_color_manual(values = c("cluster1" = "#F56E47", "cluster2" = "#97BF04", "cluster3" = "#772B59"),
#                      name = "cluster", labels = c("1", "2", "3")) +
#   scale_fill_manual(values = c("cluster1" = "#F56E47", "cluster2" = "#97BF04", "cluster3" = "#772B59"),
#                     name = "cluster", labels = c("1", "2", "3")) +
#   theme_bw() + 
#   coord_cartesian() +
#   xlim(-5, 15) +
#   ylim(-1.8, 2.5) +
#   theme(panel.border = element_rect(color="black"), 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         axis.title = element_text(size=15),
#         axis.text = element_text(size=15),
#         legend.text = element_text(size=15),
#         legend.title = element_text(size=15),
#         legend.position = c(.9, .8))

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3g.genes.pdf", 
    width=5, height=4)
scatterplot3d(x = mca_genes_df$X,
                    y = -mca_genes_df$Z,
                    z = -mca_genes_df$Y, pch = 16, 
                    color=palette2[mca_genes_df$cluster],
                    angle = 40, grid = F, box = T,
                    scale.y = 0.5)
dev.off()


# # 7. PC2 vs. PC3
# 
# ## 7.1. 3d plot of the variables (marks)
# p3 <- ggplot(mca_vars_df, aes(x=Y, y=Z, fill=group, label=mark)) +
#   geom_point(aes(size=frequency), color="black",shape=21) +
#   scale_color_manual(values = palette) +
#   scale_fill_manual(values = palette) +
#   geom_text_repel(aes(alpha=alpha, color=group), size=4.5) +
#   theme_bw() +
#   coord_cartesian() +
#   xlim(-2.5, 2.5) +
#   ylim(-3, 4.5) +
#   guides(alpha=F, color=F, fill=F) +
#   theme(panel.border = element_rect(color="black"), 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         axis.title = element_text(size=15),
#         axis.text = element_text(size=15),
#         legend.text = element_text(size=15),
#         legend.title = element_text(size=15),
#         legend.position = "top")
# 
# # p3 <- ggMarginal(p3,
# #                  type = "histogram",
# #                  groupFill = F,
# #                  fill="white",
# #                  color="white")
# 
# ## 7.2. 3d plot of the genes
# p4 <- ggplot(mca_genes_df, aes(x=Y, y=Z, color=cluster, fill=cluster)) +
#   #geom_point(size=3) +
#   stat_density_2d() +
#   scale_color_manual(values = c("cluster1" = "#F56E47", "cluster2" = "#97BF04", "cluster3" = "#772B59"),
#                      name = "cluster", labels = c("1", "2", "3")) +
#   scale_fill_manual(values = c("cluster1" = "#F56E47", "cluster2" = "#97BF04", "cluster3" = "#772B59"),
#                     name = "cluster", labels = c("1", "2", "3")) +
#   theme_bw() + 
#   coord_cartesian() +
#   xlim(-2.5, 2.5) +
#   ylim(-3, 4.5) +
#   theme(panel.border = element_rect(color="black"), 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         axis.title = element_text(size=15),
#         axis.text = element_text(size=15),
#         legend.text = element_text(size=15),
#         legend.title = element_text(size=15),
#         legend.position = "top")
# 
# # p4 <- ggMarginal(p4,
# #                  type = "histogram",
# #                  groupColour = F,
# #                  groupFill = T,
# #                  color="black")
# 
# 
# main.p2 <- plot_grid(plotlist = list(p3, p4), nrow=1, ncol=2, rel_heights = c(1, 1))
# 
# 
# pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3g2.pdf", 
#     width=16, height=4.5)
# main.p2
# dev.off()
# 
# 
# 
# # 8. retrieve clusters
# # cl1 <- as.data.frame(rownames(mca_genes_df[mca_genes_df$cluster==1, ]))
# # cl2 <- as.data.frame(rownames(mca_genes_df[mca_genes_df$cluster==2, ]))
# # cl3 <- as.data.frame(rownames(mca_genes_df[mca_genes_df$cluster==3, ]))
# 
# # write.table(cl1, "~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.1.txt", row.names = F, col.names = F, quote=F)
# # write.table(cl2, "~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.2.txt", row.names = F, col.names = F, quote=F)
# # write.table(cl3, "~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.3.txt", row.names = F, col.names = F, quote=F)
# 
# 
# 
# 
