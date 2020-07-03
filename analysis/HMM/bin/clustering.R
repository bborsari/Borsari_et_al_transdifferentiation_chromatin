.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")



#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option( c("-n", "--n_states"), type = "numeric",
               help = "Number of states." ),
  
  make_option( c("-d", "--distance_method"), 
               help = "Distance method used by TraMineR.")
  
  # make_option( c("-m", "--clustering_method"), 
  #              help = "Clustering method.")
  
  
)


parser <- OptionParser(
  usage = "%prog [options] files", 
  option_list=option_list,
  description = "\nWrapper for HMM"
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options


#************
# LIBRARIES *
#************

library(RColorBrewer)
library(cowplot)
library(reshape2)
library(dplyr)
library(tidyr)
library(TraMineR)
library(ggplotify)
library(pheatmap)


#********
# BEGIN *
#********

# 0. set working directory
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/")


# 1. read metadata
metadata <- read.table("H3K4me3/QN.merged/metadata.tsv", h=T, 
                       sep = "\t", stringsAsFactors = F)

# 2. read expression matrix of 12448 genes
expression.matrix <- read.table("expression/QN.merged/selected.genes.rep.2.3.after.QN.merged.tsv",
                                h=T, sep="\t")

# 3. retrieve set of not expressed genes
not.expressed.genes <- read.table("expression/silent.genes.txt", h=F, sep="\t", 
                                  stringsAsFactors = F)
not.expressed.genes <- not.expressed.genes$V1


# 4. retrieve set of DE genes
DE.genes <- read.table("expression/QN.merged/expression.matrix.tsv", h=T, sep="\t")
DE.genes <- rownames(DE.genes)


# 5. retrieve set of stably expressed genes 
stably.expressed.genes <- setdiff(rownames(expression.matrix), 
                                  c(not.expressed.genes, DE.genes))


# 6. change wd
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/HMM/")
i=opt$n_states

# 7.1. read emission matrix
m1 <- read.table(paste0("marks/HMM.", i, ".response.matrix.tsv"), 
                 h=T, sep="\t")

# 7.2. discard sd values
m1 <- m1[, paste0("Re", 1:9, "..Intercept.")]

# 7.3. update colnames
colnames(m1) <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3",
                  "H3K36me3", "H4K20me1", "H3K9me3", "H3K27me3")

# 7.4. sort states according to mean row value
m1$mean <- apply(m1, 1, mean)
m1 <- m1[order(m1$mean), ]
m1$mean <- NULL
z <- paste0(letters[1:i], ":", gsub("St", "", rownames(m1)))
names(z) <- rownames(m1)


# 8. read gene state sequences
m3 <- read.table(paste0("marks/HMM.", i, ".gene.matrix.tsv"), 
                 h=T, sep="\t")
m4 <- as.data.frame(t(apply(m3, 1, function(x){z[paste0("St",x)]})))
colnames(m4) <- colnames(m3)


# 9. create seqdef object
m3.seq <- seqdef(m4)


# 10. compute distance object
if (opt$distance_method %in% c("OM", "LCS", "LCP", "RLCP")) {
  
  m3.OM <- as.data.frame(seqdist(m3.seq, method = opt$distance_method,
                                 sm = "TRATE"))
  
} else {
  
  m3.OM <- as.data.frame(seqdist(m3.seq, method = opt$distance_method))
  
  
}


colnames(m3.OM) <- rownames(m3)
rownames(m3.OM) <- colnames(m3.OM)


# 11. remove numeric values from states
m4 <- as.data.frame(apply(m4, 2, function(x){gsub("\\:.*","", x)}))


# 12. convert letter states to numbers
m5 <- as.data.frame(apply(m4, 2, function(x){x <- as.numeric(as.factor(x))}))
rownames(m5) <- rownames(m4)


# 13. prepare palette
palette <- rev(rainbow(i))


# 14. compute number of states per gene
m5$n_states <- apply(m5, 1, function(x){length(table(x))})


# 15. make plots
lop <- list()
k=1

clustering_methods <- c("ward.D", "ward.D2", "single", 
                        "complete", "average", "mcquitty", 
                        "median","centroid")

for ( cm in clustering_methods ) {
  
  lop[[k]] <- as.grob(pheatmap(m5[rownames(m5) %in% not.expressed.genes, 1:12],
                               cluster_cols = F, color = palette,
                               clustering_distance_rows = as.dist(m3.OM[rownames(m3.OM) %in% not.expressed.genes,
                                                                        colnames(m3.OM) %in% not.expressed.genes]),
                               show_rownames = F,
                               border_color = NA,
                               legend = F,
                               clustering_method = cm,
                               main = paste0("not exp. genes - ", cm))$gtable)
  
  k <- k+1
  
  
}


for ( cm in clustering_methods ) {
  
  lop[[k]] <- as.grob(pheatmap(m5[rownames(m5) %in% stably.expressed.genes, 1:12],
                               cluster_cols = F, color = palette,
                               clustering_distance_rows = as.dist(m3.OM[rownames(m3.OM) %in% stably.expressed.genes,
                                                                        colnames(m3.OM) %in% stably.expressed.genes]),
                               show_rownames = F,
                               border_color = NA,
                               legend = F,
                               clustering_method = cm,
                               main = paste0("stably exp. genes - ", cm))$gtable)
  
  k <- k+1
  
  
}


for ( cm in clustering_methods ) {
  
  lop[[k]] <- as.grob(pheatmap(m5[rownames(m5) %in% DE.genes, 1:12],
                               cluster_cols = F, color = palette,
                               clustering_distance_rows = as.dist(m3.OM[rownames(m3.OM) %in% DE.genes,
                                                                        colnames(m3.OM) %in% DE.genes]),
                               show_rownames = F,
                               border_color = NA,
                               legend = F,
                               clustering_method = cm,
                               main = paste0("DE genes - ", cm))$gtable)
  
  k <- k+1
  
  
}



# 16. save plots
pdf(paste0("~/public_html/Borsari_et_al_transdifferentiation_chromatin/HMM/clustering/HMM.",
           opt$n_states, ".", opt$distance_method, ".pdf"),
    width = 22, height = 15)
plot_grid(plotlist = lop, ncol = 8)
dev.off()




# # 15.1. DE & variable genes
# for (class in c("upregulation", "downregulation", "peaking", "bending")) {
#   
#   x <- rownames(metadata[metadata$final_class == class, ])
#   tmp <- m5[rownames(m5) %in% x & m5$n_states > 1, 1:12]
#   
#   lop[[k]] <- as.grob(pheatmap(tmp, cluster_cols = F, color = palette,
#                                clustering_distance_rows = as.dist(m3.OM[rownames(m3.OM) %in% rownames(tmp),
#                                                                 colnames(m3.OM) %in% rownames(tmp)]), 
#                                show_rownames = F,
#                                clustering_method = opt$clustering_method,
#                                border_color = NA,
#                                legend = F,
#                                main = paste0(class, " - variable genes"))$gtable)
#   
#   k <- k+1
#   
# }
# 
# 
# # 15.2. DE & stable genes
# for (class in c("upregulation", "downregulation", "peaking", "bending")) {
#   
#   x <- rownames(metadata[metadata$final_class == class, ])
#   tmp <- m5[rownames(m5) %in% x & m5$n_states == 1, 1:12]
#   
#   lop[[k]] <- as.grob(pheatmap(tmp, cluster_cols = F, color = palette,
#                                clustering_distance_rows = as.dist(m3.OM[rownames(m3.OM) %in% rownames(tmp),
#                                                                         colnames(m3.OM) %in% rownames(tmp)]), 
#                                show_rownames = F,
#                                clustering_method = opt$clustering_method,
#                                border_color = NA,
#                                legend = F,
#                                main = paste0(class, " - stable genes"))$gtable)
#   
#   k <- k+1
#   
# }
# 
# 
# # 15.3. not_expressed genes
# lop[[k]] <- as.grob(pheatmap(m5[rownames(m5) %in% not.expressed.genes, 1:12], 
#                              cluster_cols = F, color = palette,
#                              clustering_distance_rows = as.dist(m3.OM[rownames(m3.OM) %in% not.expressed.genes,
#                                                                       colnames(m3.OM) %in% not.expressed.genes]), 
#                              show_rownames = F,
#                              border_color = NA,
#                              legend = F,
#                              clustering_method = opt$clustering_method,
#                              main = "not expressed genes")$gtable)
# k <- k+1
# 
# 
# # 15.4. stably expressed genes
# lop[[k]] <- as.grob(pheatmap(m5[rownames(m5) %in% stably.expressed.genes, 1:12], 
#                              cluster_cols = F, color = palette,
#                              clustering_distance_rows = as.dist(m3.OM[rownames(m3.OM) %in% stably.expressed.genes,
#                                                                       colnames(m3.OM) %in% stably.expressed.genes]), 
#                              show_rownames = F,
#                              border_color = NA,
#                              legend = F,
#                              clustering_method = opt$clustering_method,
#                              main = "stably expressed genes")$gtable)
# k <- k+1
# 
# 
# # 15.4. stably expressed genes
# lop[[k]] <- as.grob(pheatmap(m5[rownames(m5) %in% DE.genes, 1:12], 
#                              cluster_cols = F, color = palette,
#                              clustering_distance_rows = as.dist(m3.OM[rownames(m3.OM) %in% DE.genes,
#                                                                       colnames(m3.OM) %in% DE.genes]), 
#                              show_rownames = F,
#                              border_color = NA,
#                              legend = F,
#                              clustering_method = opt$clustering_method,
#                              main = "DE genes")$gtable)
# 
# 
# 
