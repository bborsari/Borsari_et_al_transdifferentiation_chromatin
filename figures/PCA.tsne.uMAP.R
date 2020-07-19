.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")



#************
# LIBRARIES *
#************

library(plotly)
library(Rtsne)
library(umap)


# 1. set wd
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/")


# 2. the marks we're analyzing
marks <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3", "H3K36me3",
           "H4K20me1", "H3K9me3", "H3K27me3")

# 3. read expression matrix
m <- read.table("expression/QN.merged/expression.matrix.tsv", h=T, sep="\t")
colnames(m) <- paste0(colnames(m), "_expression")


# 4. read marks' matrices
for ( x in marks ) {
  
  tmp <- read.table(paste0(x, "/QN.merged/", x, ".matrix.tsv"), h=T, sep="\t")
  stopifnot(identical(rownames(m), rownames(tmp)))
  colnames(tmp) <- paste0(colnames(tmp), "_", x)
  
  m <- cbind(m, tmp)
  
  
}


# clusters
clusters <- data.frame(stringsAsFactors = F)
for (k in 1:3) {
  
  tmp <- read.table(paste0("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.", k, ".txt"),
                    h=F, stringsAsFactors = F, sep="\t")
  tmp$cluster <- paste0("cluster", k)
  clusters <- rbind(clusters, tmp)
  
}
colnames(clusters)[1] <- "gene_id"
rownames(clusters) <- clusters$gene_id
clusters <- clusters[rownames(m), ]

palette <- c("cluster1" = "#F56E47", "cluster2" = "#97BF04", "cluster3" = "#772B59")


# 5. make PCA with expression and marks
pca.obj1 <- prcomp(m, center = T, scale = F)
summary(pca.obj1)
pca.m1 <- as.data.frame(pca.obj1$x[, 1:3])

pca.1 <- plot_ly(pca.m1, x = ~PC1, y = ~PC3, z = ~PC2, size = 1, 
                 color = as.factor(clusters$cluster), 
                 colors = palette) %>%
  add_markers()
htmlwidgets::saveWidget(as_widget(pca.1),
                        "~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/pca.expression.marks.html")


# 6. make PCA with only marks
pca.obj2 <- prcomp(m[, 13:120], center = T, scale = F)
summary(pca.obj2)
pca.m2 <- as.data.frame(pca.obj2$x[, 1:3])

pca.2 <- plot_ly(pca.m2, x = ~PC1, y = ~PC3, z = ~PC2, size = 1,
                 color =  as.factor(clusters$cluster),
                 colors = palette) %>%
  add_markers()
htmlwidgets::saveWidget(as_widget(pca.2),
                        "~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/pca.marks.html")


# 7. make tsne with expression and marks
set.seed(123)
tsne.obj1 <- Rtsne(m, dims = 3, perplexity=30)
tsne.m1 <- as.data.frame(tsne.obj1$Y)
rownames(tsne.m1) <- rownames(m)
tsne.1 <- plot_ly(tsne.m1, x = ~V1, y=~V3, z=~V2, size = 1,
                  color = as.factor(clusters$cluster), 
                  colors = palette) %>%
  add_markers()
htmlwidgets::saveWidget(as_widget(tsne.1),
                        "~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/tsne.expression.marks.html")


# 8. make tsne with only marks
set.seed(123) 
tsne.obj2 <- Rtsne(m[, 13:120], dims = 3, perplexity=30)
tsne.m2 <- as.data.frame(tsne.obj2$Y)
rownames(tsne.m2) <- rownames(m)
tsne.2 <- plot_ly(tsne.m2, x = ~V1, y=~V3, z=~V2, size = 1,
                  color = as.factor(clusters$cluster),
                  colors = palette) %>%
  add_markers()
htmlwidgets::saveWidget(as_widget(tsne.2),
                        "~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/tsne.marks.html")


# 9. make umap with expression and marks
set.seed(123)
umap.obj1 <- umap(m)
umap.m1 <- as.data.frame(umap.obj1$layout)
umap.1 <- plot_ly(umap.m1, x = ~V1, y=~V2, size = 1,
                  color = as.factor(clusters$cluster),
                  colors = palette) %>%
  add_markers()
htmlwidgets::saveWidget(as_widget(umap.1),
                        "~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/umap.expression.marks.html")



# 10. make umap with marks only
set.seed(123)
umap.obj2 <- umap(m[, 13:120])
umap.m2 <- as.data.frame(umap.obj2$layout)
umap.2 <- plot_ly(umap.m2, x = ~V1, y=~V2, size = 1,
                  color = as.factor(clusters$cluster),
                  colors = palette) %>%
  add_markers()
htmlwidgets::saveWidget(as_widget(umap.2),
                        "~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/umap.marks.html")

