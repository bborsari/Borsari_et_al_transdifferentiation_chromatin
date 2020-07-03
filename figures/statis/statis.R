.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

#************
# LIBRARIES *
#************

library(ade4)
library(ggplot2)

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/")

#********
# BEGIN *
#********

# 1. read dataframes

expression.matrix <- read.table("expression/QN.merged/selected.genes.rep.2.3.after.QN.merged.tsv",
                                h=T, sep="\t")
H3K27ac.matrix <- read.table("H3K27ac/QN.merged/H3K27ac.matrix.after.QN.merged.tsv",
                             h=T, sep="\t")
H3K9ac.matrix <- read.table("H3K9ac/QN.merged/H3K9ac.matrix.after.QN.merged.tsv",
                            h=T, sep="\t")
H4K20me1.matrix <- read.table("H4K20me1/QN.merged/H4K20me1.matrix.after.QN.merged.tsv",
                              h=T, sep="\t")
H3K4me3.matrix <- read.table("H3K4me3/QN.merged/H3K4me3.matrix.after.QN.merged.tsv",
                             h=T, sep="\t")
H3K4me1.matrix <- read.table("H3K4me1/QN.merged/H3K4me1.matrix.after.QN.merged.tsv",
                             h=T, sep="\t")
H3K36me3.matrix <- read.table("H3K36me3/QN.merged/H3K36me3.matrix.after.QN.merged.tsv",
                              h=T, sep="\t")
H3K4me2.matrix <- read.table("H3K4me2/QN.merged/H3K4me2.matrix.after.QN.merged.tsv",
                             h=T, sep="\t")
H3K9me3.matrix <- read.table("H3K9me3/QN.merged/H3K9me3.matrix.after.QN.merged.tsv",
                             h=T, sep="\t")
H3K27me3.matrix <- read.table("H3K27me3/QN.merged/H3K27me3.matrix.after.QN.merged.tsv",
                              h=T, sep="\t")

# 2. check order of rownames
expression.matrix <- expression.matrix[rownames(H3K27ac.matrix), ]

stopifnot(identical(rownames(expression.matrix),
                    rownames(H3K27ac.matrix)))
stopifnot(identical(rownames(expression.matrix),
                    rownames(H3K9ac.matrix)))
stopifnot(identical(rownames(expression.matrix),
                    rownames(H4K20me1.matrix)))
stopifnot(identical(rownames(expression.matrix),
                    rownames(H3K4me3.matrix)))
stopifnot(identical(rownames(expression.matrix),
                    rownames(H3K4me1.matrix)))
stopifnot(identical(rownames(expression.matrix),
                    rownames(H3K36me3.matrix)))
stopifnot(identical(rownames(expression.matrix),
                    rownames(H3K4me2.matrix)))
stopifnot(identical(rownames(expression.matrix),
                    rownames(H3K9me3.matrix)))
stopifnot(identical(rownames(expression.matrix),
                    rownames(H3K27me3.matrix)))

# 3. get transpose & scaled matrices
expression.matrix <- as.data.frame(t(expression.matrix))
H3K27ac.matrix <- as.data.frame(t(H3K27ac.matrix))
H3K9ac.matrix <- as.data.frame(t(H3K9ac.matrix))
H4K20me1.matrix <- as.data.frame(t(H4K20me1.matrix))
H3K4me3.matrix <- as.data.frame(t(H3K4me3.matrix))
H3K4me1.matrix <- as.data.frame(t(H3K4me1.matrix))
H3K36me3.matrix <- as.data.frame(t(H3K36me3.matrix))
H3K4me2.matrix <- as.data.frame(t(H3K4me2.matrix))
H3K9me3.matrix <- as.data.frame(t(H3K9me3.matrix))
H3K27me3.matrix <- as.data.frame(t(H3K27me3.matrix))




###----
### STATIS
###----


# test <- rbind(expression.matrix,
#               H3K27ac.matrix,
#               H3K9ac.matrix,
#               H4K20me1.matrix,
#               H3K4me3.matrix,
#               H3K4me1.matrix,
#               H3K36me3.matrix,
#               H3K4me2.matrix,
#               H3K9me3.matrix,
#               H3K27me3.matrix)

test <- rbind(scale(expression.matrix),
              scale(H3K27ac.matrix))
test2 <- test[, 1:100]


# lab <- as.factor(rep(rownames(expression.matrix), 10))

lab <- as.factor(rep(rownames(expression.matrix), 2))


# info <- rep(c("expression",
#               "H3K27ac",
#               "H3K9ac",
#               "H4K20me1",
#               "H3K4me3",
#               "H3K4me1",
#               "H3K36me3",
#               "H3K4me2",
#               "H3K9me3",
#               "H3K27me3"), each = 12)

kta1 <- ktab.within(withinpca(test2, lab, scann = FALSE, scaling = "partial"))
kta2 <- ktab.within(withinpca(test2, lab, scann = FALSE, scaling = "total"))

statis1 <- statis(kta1, scann = FALSE)
statis2 <- statis(kta2, scann = FALSE)


plot(statis1)
plot(statis2)


df1 <- statis1$C.li
ggplot(df1, aes(C2, C3)) + geom_point()


df2 <- statis2$C.Co
ggplot(df2, aes(C2, C3)) + geom_point()

save.image(file="/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/figures/statis.expression.H3K27ac.RData")
