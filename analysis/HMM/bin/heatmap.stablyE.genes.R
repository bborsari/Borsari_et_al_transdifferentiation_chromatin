.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")



#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option( c("-n", "--n_states"), type = "numeric",
               help = "Number of states." )
  
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


#************
# FUNCTIONS *
#************


mypdf = function(pdfname, mypattern = "MYTEMPPNG", ...) {
  fname = paste0(mypattern, "%05d.png")
  gpat = paste0(mypattern, ".*\\.png")
  takeout = list.files(path = tempdir(), pattern = gpat, full.names = TRUE)
  if (length(takeout) > 0)
    file.remove(takeout)
  pngname = file.path(tempdir(), fname)
  png(pngname, ...)
  return(list(pdfname = pdfname, mypattern = mypattern))
}


mydev.off = function(pdfname, mypattern, copts = "") {
  dev.off()
  gpat = paste0(mypattern, ".*\\.png")
  pngs = list.files(path = tempdir(), pattern = gpat, full.names = TRUE)
  mystr = paste(pngs, collapse = " ", sep = "")
  system(sprintf("convert %s -quality 100 %s %s", mystr, pdfname, copts))
}



#********
# BEGIN *
#********

# 1. set working directory
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/")


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


# 6. create a df with gene label (not expressed, stably expressed, DE)
x <- data.frame(gene_id = c(DE.genes, stably.expressed.genes, not.expressed.genes),
                group = c(rep("DE", length(DE.genes)), 
                          rep("stably e.", length(stably.expressed.genes)),
                          rep("not e.", length(not.expressed.genes))))

# 7. add clusters info for DE genes
clusters <- data.frame(stringsAsFactors = F)
for (k in 1:3) {
  
  tmp <- read.table(paste0("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.", k, ".txt"),
                    h=F, stringsAsFactors = F, sep="\t")
  tmp$cluster <- paste0("cluster", k)
  clusters <- rbind(clusters, tmp)
  
}
colnames(clusters)[1] <- "gene_id"

x <- merge(x, clusters, all.x = T)


# 8. add class (upreg/downreg/peaking/bending) info for DE genes
metadata <- read.table("H3K4me3/QN.merged/metadata.tsv", h=T, sep="\t", stringsAsFactors = F)
metadata$gene_id <- rownames(metadata)

x <- merge(x, metadata[, c("gene_id", "final_class")], by = "gene_id", all.x = T)
# x$cluster <- ifelse(is.na(x$cluster), as.character(x$group), as.character(x$cluster))
# x$final_class <- ifelse(is.na(x$final_class), as.character(x$group), as.character(x$final_class))


# 9. set HMM wd
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/HMM/")
i=opt$n_states

# 10. read emission matrix
m1 <- read.table(paste0("marks/HMM.", i, ".response.matrix.tsv"), 
                 h=T, sep="\t")

# 11. discard sd values
m1 <- m1[, paste0("Re", 1:9, "..Intercept.")]

# 12. update colnames
colnames(m1) <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3",
                  "H3K36me3", "H4K20me1", "H3K9me3", "H3K27me3")

# 13. sort states according to mean row value
m1$mean <- apply(m1, 1, mean)
m1 <- m1[order(m1$mean), ]
m1$mean <- NULL
z <- paste0(letters[1:i], ":", gsub("St", "", rownames(m1)))
names(z) <- rownames(m1)


# 14. read gene state sequences
m3 <- read.table(paste0("marks/HMM.", i, ".gene.matrix.tsv"), 
                 h=T, sep="\t")
m4 <- as.data.frame(t(apply(m3, 1, function(x){z[paste0("St",x)]})))
colnames(m4) <- colnames(m3)


# 15. create seqdef object
m3.seq <- seqdef(m4)


# 16. compute distance object
m3.OM <- as.data.frame(seqdist(m3.seq, method = "NMS"))

colnames(m3.OM) <- rownames(m3)
rownames(m3.OM) <- colnames(m3.OM)


# 17. remove numeric values from states
m4 <- as.data.frame(apply(m4, 2, function(x){gsub("\\:.*","", x)}))


# 18. convert letter states to numbers
m5 <- as.data.frame(apply(m4, 2, function(x){x <- as.numeric(as.factor(x))}))
rownames(m5) <- rownames(m4)


# 19. prepare palette
palette <- rev(rainbow(i))


# 20. compute number of states per gene
m5$n_states <- apply(m5, 1, function(x){length(table(x))})


# 21. make plots
rownames(x) <- x$gene_id
stopifnot(identical(rownames(x), rownames(m5)))
x$cluster <- factor(x$cluster)



clustering_methods <- c("ward.D", "ward.D2", "single", 
                        "complete", "average", "mcquitty", 
                        "median","centroid")

for (cm in clustering_methods) {
  
  res <- mypdf(pdfname = paste0("~/public_html/Borsari_et_al_transdifferentiation_chromatin/HMM/n_states/",
                                opt$n_states, "/heatmap.", cm, ".stablyE.genes.pdf"), 
               res = 300, height = 17, width = 14, units = "cm")
  
  pheatmap(m5[rownames(m5) %in% stably.expressed.genes, 1:12],
           cluster_cols = F, color = palette,
           clustering_distance_rows = as.dist(m3.OM[rownames(m3.OM) %in% stably.expressed.genes,
                                                    colnames(m3.OM) %in% stably.expressed.genes]),
           annotation_row = x[rownames(x) %in% stably.expressed.genes, 2:4],
           show_rownames = F,
           border_color = NA,
           legend = F,
           clustering_method = cm,
           annotation_colors = list(group = c("DE" = "#f768a1", 
                                              "stably e." = "#bdbdbd", 
                                              "not e." = "#737373"),
                                    cluster = c("cluster1" = "#F56E47", 
                                                "cluster2" = "#97BF04", 
                                                "cluster3" = "#772B59"),
                                    final_class = c("downregulation" = "#2d7f89", 
                                                    "bending" = "#7acbd5",
                                                    "upregulation" = "#89372d",
                                                    "peaking" = "#d5847a")),
           main = cm)
  
  mydev.off(pdfname = res$pdfname, mypattern = res$mypattern)
  
  
  
}

