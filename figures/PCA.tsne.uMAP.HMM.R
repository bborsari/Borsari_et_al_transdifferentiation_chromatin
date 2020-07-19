.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")



#************
# LIBRARIES *
#************

library(plotly)
library(Rtsne)
library(umap)


# 1. set wd
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/HMM/marks")


# 2. read HMM matrix
m <- read.table("HMM.5.gene.matrix.tsv", h=T, sep="\t")


# 3. read clusters df
clusters <- data.frame(stringsAsFactors = F)
for (k in 1:3) {
  
  tmp <- read.table(paste0("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.", k, ".txt"),
                    h=F, stringsAsFactors = F, sep="\t")
  tmp$cluster <- paste0("cluster", k)
  clusters <- rbind(clusters, tmp)
  
}
colnames(clusters)[1] <- "gene_id"
rownames(clusters) <- clusters$gene_id


# 4. keep only DE genes
m <- m[rownames(m) %in% rownames(clusters), ]
m <- m[rownames(clusters), ]
stopifnot(identical(rownames(m), rownames(clusters)))

palette <- c("cluster1" = "#F56E47", "cluster2" = "#97BF04", "cluster3" = "#772B59")


# 5. make PCA with HMM matrix
pca.obj1 <- prcomp(m, center = T, scale = F)
summary(pca.obj1)
pca.m1 <- as.data.frame(pca.obj1$x[, 1:3])

pca.1 <- plot_ly(pca.m1, x = ~PC1, y = ~PC3, z = ~PC2, size = 1, 
                 color = as.factor(clusters$cluster), 
                 colors = palette) %>%
  add_markers()
htmlwidgets::saveWidget(as_widget(pca.1),
                        "~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/pca.HMM.marks.html")



# 6. make tsne with HMM matrix
set.seed(123)
tsne.obj1 <- Rtsne(unique(m), dims = 3, perplexity=30)
tsne.m1 <- as.data.frame(tsne.obj1$Y)
rownames(tsne.m1) <- rownames(unique(m))
tsne.1 <- plot_ly(tsne.m1, x = ~V1, y=~V3, z=~V2, size = 0.5,
                  color = as.factor(clusters[rownames(clusters) %in% rownames(unique(m)), "cluster"]), 
                  colors = palette) %>%
  add_markers()
htmlwidgets::saveWidget(as_widget(tsne.1),
                        "~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/tsne.HMM.marks.html")



# 7. make umap with HMM matrix
set.seed(123)
umap.obj1 <- umap(m)
umap.m1 <- as.data.frame(umap.obj1$layout)
umap.1 <- plot_ly(umap.m1, x = ~V1, y=~V2, size = 0.5,
                  color = as.factor(clusters$cluster),
                  colors = palette) %>%
  add_markers()
htmlwidgets::saveWidget(as_widget(umap.1),
                        "~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/umap.HMM.marks.html")



