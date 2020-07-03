.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


#************
# LIBRARIES *
#************

library(gridExtra)
library(xtable)




#********
# BEGIN *
#********


# 1. set working directory
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


# 2. the marks we're analyzing
marks <- c("H3K27ac", "H3K9ac", "H4K20me1", "H3K4me3", "H3K4me1", "H3K36me3",
           "H3K4me2", "H3K9me3", "H3K27me3")


# 3. read expression matrix of 12448 genes
expression.matrix <- read.table("expression/QN.merged/selected.genes.rep.2.3.after.QN.merged.tsv",
                                h=T, sep="\t")


# 4. retrieve set of not expressed genes
not.expressed.genes <- read.table("expression/silent.genes.txt", h=F, sep="\t", 
                                  stringsAsFactors = F)
not.expressed.genes <- not.expressed.genes$V1


# 5. retrieve set of DE genes
DE.genes <- read.table("expression/QN.merged/expression.matrix.tsv", h=T, sep="\t")
DE.genes <- rownames(DE.genes)


# 6. retrieve set of stably expressed genes
stably.expressed.genes <- setdiff(rownames(expression.matrix),
                                  c(not.expressed.genes, DE.genes))


# 7. read dfs of marked genes and store:
v1 <- c() # total number of marked genes
v2 <- c() # DE & marked
v3 <- c() # stably.expressed & marked
v4 <- c() # not.expressed & marked

for ( i in 1:9 ) {
  
  x <- read.table(paste0(marks[i], "/QN.merged/all.genes.intersecting.peaks.tsv"),
                  h=F, sep="\t", stringsAsFactors = F)
  
  # 7.1. # total number of marked genes
  v1 <- c(v1, nrow(x))
  
  # 7.2. DE & marked
  v2 <- c(v2, length(intersect(x$V1, DE.genes)))
  
  # 7.3. stably.expressed & marked
  v3 <- c(v3, length(intersect(x$V1, stably.expressed.genes)))
  
  # 7.4. not.expressed & marked
  v4 <- c(v4, length(intersect(x$V1, not.expressed.genes)))
  
}


# 8. prepare df for table
df <- data.frame(marks = marks,
                 marked.genes = v1,
                 DE.marked = v2,
                 stably.expressed.marked = v3,
                 not.expressed.marked = v4)


# 9. reorder marks by the total number of marked genes
df <- df[order(df$marked.genes, decreasing = T), ]
df$marked.genes <- NULL


# 10. group DE and stably.expressed = expressed
df$expressed.marked <- df$DE.marked + df$stably.expressed.marked
df$DE.marked <- NULL
df$stably.expressed.marked <- NULL


# 11. obtain expressed and not_marked
df$expressed.not.marked <- 10696 - df$expressed.marked
df$expressed.not.marked <- paste0(df$expressed.not.marked, " (", 
                                  round(((df$expressed.not.marked / 10696)*100), 2), "%)")
df$expressed.marked <- paste0(df$expressed.marked, " (", 
                              round(((df$expressed.marked / 10696)*100), 2), "%)")


# 12. obtain not_expressed and not_marked
df$not.expressed.not.marked <- 1552 - df$not.expressed.marked
df$not.expressed.not.marked <- paste0(df$not.expressed.not.marked, " (", 
                                  round(((df$not.expressed.not.marked / 1552)*100), 2), "%)")
df$not.expressed.marked <- paste0(df$not.expressed.marked, " (", 
                                  round(((df$not.expressed.marked / 1552)*100), 2), "%)")



# 13. reorder cols df
df <- df[, c("marks", "not.expressed.not.marked", "not.expressed.marked",
             "expressed.not.marked", "expressed.marked")]


# 14. table in latex format
print(xtable(df), include.rownames=FALSE)
