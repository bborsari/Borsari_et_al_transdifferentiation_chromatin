.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(reshape2)
library(ggalluvial)
library(viridis)


#************
# FUNCTIONS *
#************

function1 <- function(mark) {
  
  df <- read.table(paste0(mark, "/QN.merged/ant.del.analysis/increases/", 
                          "expression.", mark, ".dynamics.tsv"),
                   h=T, sep="\t")
  
  df <- df[, c("group_25", "group_50", "group_75", "group_100")]
  
  df$gene_id <- rownames(df)
  
  df.melt <- melt(df, id.vars = "gene_id")
  df.melt$value <- gsub("concomitant", "co-occurrent", df.melt$value)
  df.melt$value <- factor(df.melt$value, levels = c("anticipated", "co-occurrent", "delayed"))

  df.melt$mark <- mark
  
  df.melt <- df.melt[df.melt$variable == "group_25", c("gene_id", "value")]
  colnames(df.melt)[2] <- mark
  
  return(df.melt)
  
  
}

#********
# BEGIN *
#********

# 1. set working directory
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


# 2. define palette
palette <- c("anticipated" = "#bdbdbd",
             "co-occurrent" = "#737373",
             "delayed" = "#403734")

# 3. the marks we're analyzing
marks <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3", "H3K36me3", "H4K20me1")


# 4. make unique df
x <- function1(mark=marks[1])
for ( i in 2:7 ){
  
  x <- merge(x, function1(mark=marks[i]), by = "gene_id", all=T)
  
}


# 5. convert categorical data to numerical data
z <- x$gene_id
x <- as.matrix(x)
x[x=="anticipated"] <- 1
x[x=="co-occurrent"] <- 2
x[x=="delayed"] <- 3
x <- as.data.frame(x)
rownames(x) <- x$gene_id
x$gene_id <- NULL


# 6. Remove rows with more than 1 NA
nna.th <- 7

sel <- c()
for (i in 1:nrow(x)){
  ss <- t(x[i, ])
  if(sum(!is.na(ss)) >= nna.th) {
    sel <- c(sel, i)
  }
}
x <- x[sel,]
x <- apply(x, 2, as.numeric)
rownames(x) <- z[sel]
#tx <- t(x)
# dmat2 <- cor(tx, use = "pairwise.complete.obs")

# 7. Build distance matrix using complete cases
dmat <- c()
for (i in 1:(nrow(x))){
  for (j in 1:i){
    ss <- t(x[c(i, j), ])
    ss <- ss[complete.cases(ss), ]
    d <- as.numeric(dist(t(ss), "euclidean"))
    # d <- sum(ss[, 1]!= ss[, 2]) / nrow(ss)
    # d <- as.numeric(dist(t(ss), "euclidean"))/sqrt(nrow(ss))
    # d <- 1-cor(ss[,1], ss[,2] )
    dmat <- rbind(dmat, data.frame(i = rownames(x)[i], j = rownames(x)[j], d))
  }
}




dmat <- dcast(dmat, i ~ j, value.var = "d")
rownames(dmat) <- dmat$i
dmat$i <- NULL

tdmat <- t(dmat)
dmat[upper.tri(dmat, diag = F)] <- t(tdmat[upper.tri(tdmat, diag = F)])
dmat <- as.dist(dmat)


save.image(file = "/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/figures/figure.9u.RData")


pheatmap(x, clustering_distance_rows = dmat,
         cluster_cols = F, cluster_rows = T,
         show_rownames = F,
         color = c("1" = "#bdbdbd",
                   "2" = "#737373",
                   "3" = "#403734"),
         clustering_method = "complete",
         border_color = NA,
         cellwidth = 10)

