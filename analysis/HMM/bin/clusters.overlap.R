.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")



#************
# LIBRARIES *
#************

library(pheatmap)
library(ggplotify)


#********
# BEGIN *
#********

cluster1 <- 
  read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/GO.analysis/cluster.1.txt",
             h=F, sep="\t", stringsAsFactors = F)
cluster1 <- cluster1$V1

cluster2 <- 
  read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/GO.analysis/cluster.2.txt",
             h=F, sep="\t", stringsAsFactors = F)
cluster2 <- cluster2$V1


cluster3 <- 
  read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/GO.analysis/cluster.3.txt",
             h=F, sep="\t", stringsAsFactors = F)
cluster3 <- cluster3$V1






function1 <- function(x, stable, variable) {
  
  y <- length(intersect(x, stable)) / length(x)
  z <- length(intersect(x, variable)) / length(x)
  
  print(c(y, z))
  
  
}


function2 <- function(i) {
  
  m <- read.table(paste0("marks/HMM.", i, ".gene.matrix.tsv"), 
                   h=T, sep="\t")
  
  m$n_states <- apply(m, 1, function(x){length(table(x))})
  m$n_states <- ifelse(m$n_states < 2, "stable", "variable")
  DE.genes <- read.table("../all.marks/expression/QN.merged/expression.matrix.tsv", h=T, sep="\t")
  DE.genes <- rownames(DE.genes)
  
  m.DE <- m[rownames(m) %in% DE.genes, ]
  stable <- rownames(m.DE[m.DE$n_states == "stable", ])
  variable <- rownames(m.DE[m.DE$n_states == "variable", ])
  
  
  df <- as.data.frame(rbind(function1(x=cluster1, stable = stable, variable = variable),
                            function1(x=cluster2, stable = stable, variable = variable),
                            function1(x=cluster3, stable = stable, variable = variable)))

  colnames(df) <- c("stable", "variable")
  
  df$clusters <- c("cluster1", "cluster2", "cluster3")

  df$n_states <- i
  
  
  return(df)
  
}


a <- data.frame(stringsAsFactors = F)


for ( j in 2:20) {
  
  a <- rbind(a, function2(i = j))
  
}


lop <- list()


lop[[1]] <- as.grob(pheatmap(t(a[a$clusters == "cluster1", 1:2]), cluster_rows = F,
                     cluster_cols = F,
                     labels_col = a[a$clusters == "cluster1", 4],
                     display_numbers = T,
                     cellwidth = 25,
                     cellheight = 25,
                     border_color = NA,
                     main = "cluster1",
                     breaks = seq(0, 1, by = 0.1),
                     color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(10))$gtable)


lop[[2]] <- as.grob(pheatmap(t(a[a$clusters == "cluster2", 1:2]), cluster_rows = F,
                             cluster_cols = F,
                             labels_col = a[a$clusters == "cluster2", 4],
                             display_numbers = T,
                             cellwidth = 25,
                             cellheight = 25,
                             border_color = NA,
                             main = "cluster2",
                             breaks = seq(0, 1, by = 0.1),
                             color = colorRampPalette(rev(brewer.pal(n = 7, 
                                                                     name = "RdYlBu")))(10))$gtable)

lop[[3]] <- as.grob(pheatmap(t(a[a$clusters == "cluster3", 1:2]), cluster_rows = F,
                             cluster_cols = F,
                             labels_col = a[a$clusters == "cluster3", 4],
                             display_numbers = T,
                             cellwidth = 25,
                             cellheight = 25,
                             border_color = NA,
                             main = "cluster3",                            
                             breaks = seq(0, 1, by = 0.1),
                             color = colorRampPalette(rev(brewer.pal(n = 7, 
                                                                     name = "RdYlBu")))(10))$gtable)


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/HMM/clusters.overlap.pdf",
    width = 10)
plot_grid(plotlist = lop, nrow=3)
dev.off()
