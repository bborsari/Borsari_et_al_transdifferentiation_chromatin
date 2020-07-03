.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")



# palette <- c("H3K9ac" = "#e7298a",
#              "H3K27ac" = "#8e0152",
#              "H3K4me3" = "#f1b6da",
#              "H3K27me3" = "#253494",
#              "H3K9me3" = "#41b6c4",
#              "H3K36me3" = "#7fbc41",
#              "H4K20me1" = "#276419",
#              "H3K4me1" = "#ffbf00",
#              "H3K4me2" = "#a67c00",
#              "expression" = "white")


palette <- c("H3K9ac" = "#af4d85",
             "H3K27ac" = "#630039",
             "H3K4me3" = "#d199b9",
             "H3K27me3" = "#1d2976",
             "H3K9me3" = "#a7add4",
             "H3K36me3" = "#7fbc41",
             "H4K20me1" = "#4c7027",
             "H3K4me1" = "#e5ab00",
             "H3K4me2" = "#a67c00",
             "expression" = "white")


palette2 <- c("H3K9ac" = "#af4d85",
              "H3K27ac" = "#630039",
              "H3K4me3" = "#d199b9",
              "H3K27me3" = "#1d2976",
              "H3K9me3" = "#a7add4",
              "H3K36me3" = "#7fbc41",
              "H4K20me1" = "#4c7027",
              "H3K4me1" = "#e5ab00",
              "H3K4me2" = "#a67c00",
              "expression" = "gray")




#************
# LIBRARIES *
#************

library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(cowplot)
library(reshape2)
library(ade4)


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
expression.matrix.scaled <- as.data.frame(scale(expression.matrix))
sum(apply(expression.matrix.scaled, 1, is.na))
expression.matrix.scaled <- as.data.frame(t(na.omit(t(expression.matrix.scaled))))

H3K27ac.matrix <- as.data.frame(t(H3K27ac.matrix))
H3K27ac.matrix.scaled <- as.data.frame(scale(H3K27ac.matrix))
sum(apply(H3K27ac.matrix.scaled, 1, is.na))
H3K27ac.matrix.scaled <- as.data.frame(t(na.omit(t(H3K27ac.matrix.scaled))))

H3K9ac.matrix <- as.data.frame(t(H3K9ac.matrix))
H3K9ac.matrix.scaled <- as.data.frame(scale(H3K9ac.matrix))
sum(apply(H3K9ac.matrix.scaled,1, is.na))
H3K9ac.matrix.scaled <- as.data.frame(t(na.omit(t(H3K9ac.matrix.scaled))))

H4K20me1.matrix <- as.data.frame(t(H4K20me1.matrix))
H4K20me1.matrix.scaled <- as.data.frame(scale(H4K20me1.matrix))
sum(apply(H4K20me1.matrix.scaled,1, is.na))
H4K20me1.matrix.scaled <- as.data.frame(t(na.omit(t(H4K20me1.matrix.scaled))))

H3K4me3.matrix <- as.data.frame(t(H3K4me3.matrix))
H3K4me3.matrix.scaled <- as.data.frame(scale(H3K4me3.matrix))
sum(apply(H3K4me3.matrix.scaled,1, is.na))
H3K4me3.matrix.scaled <- as.data.frame(t(na.omit(t(H3K4me3.matrix.scaled))))

H3K4me1.matrix <- as.data.frame(t(H3K4me1.matrix))
H3K4me1.matrix.scaled <- as.data.frame(scale(H3K4me1.matrix))
sum(apply(H3K4me1.matrix.scaled,1, is.na))
H3K4me1.matrix.scaled <- as.data.frame(t(na.omit(t(H3K4me1.matrix.scaled))))

H3K36me3.matrix <- as.data.frame(t(H3K36me3.matrix))
H3K36me3.matrix.scaled <- as.data.frame(scale(H3K36me3.matrix))
sum(apply(H3K36me3.matrix.scaled,1, is.na))
H3K36me3.matrix.scaled <- as.data.frame(t(na.omit(t(H3K36me3.matrix.scaled))))

H3K4me2.matrix <- as.data.frame(t(H3K4me2.matrix))
H3K4me2.matrix.scaled <- as.data.frame(scale(H3K4me2.matrix))
sum(apply(H3K4me2.matrix.scaled,1, is.na))
H3K4me2.matrix.scaled <- as.data.frame(t(na.omit(t(H3K4me2.matrix.scaled))))

H3K9me3.matrix <- as.data.frame(t(H3K9me3.matrix))
H3K9me3.matrix.scaled <- as.data.frame(scale(H3K9me3.matrix))
sum(apply(H3K9me3.matrix.scaled,1, is.na))
H3K9me3.matrix.scaled <- as.data.frame(t(na.omit(t(H3K9me3.matrix.scaled))))

H3K27me3.matrix <- as.data.frame(t(H3K27me3.matrix))
H3K27me3.matrix.scaled <- as.data.frame(scale(H3K27me3.matrix))
sum(apply(H3K27me3.matrix.scaled,1, is.na))
H3K27me3.matrix.scaled <- as.data.frame(t(na.omit(t(H3K27me3.matrix.scaled))))







###----
### STATIS scaling each matrix separately
###----



# 9. STATIS 

my.genes <- Reduce(intersect, list(colnames(expression.matrix.scaled),
                                   colnames(H3K27ac.matrix.scaled),
                                   colnames(H3K9me3.matrix.scaled),
                                   colnames(H4K20me1.matrix.scaled),
                                   colnames(H3K4me3.matrix.scaled),
                                   colnames(H3K4me1.matrix.scaled),
                                   colnames(H3K36me3.matrix.scaled),
                                   colnames(H3K4me2.matrix.scaled),
                                   colnames(H3K9me3.matrix.scaled),
                                   colnames(H3K27me3.matrix.scaled)))


w1 <- cbind(t(expression.matrix.scaled[, my.genes]),
            t(H3K27ac.matrix.scaled[, my.genes]),
            t(H3K9ac.matrix.scaled[, my.genes]),
            t(H4K20me1.matrix.scaled[, my.genes]),
            t(H3K4me3.matrix.scaled[, my.genes]),
            t(H3K4me1.matrix.scaled[, my.genes]),
            t(H3K36me3.matrix.scaled[, my.genes]),
            t(H3K4me2.matrix.scaled[, my.genes]),
            t(H3K9me3.matrix.scaled[, my.genes]),
            t(H3K27me3.matrix.scaled[, my.genes]))

# lab <- c("expression" = 100, "H3K27ac" = 100, "H3K9ac" = 100,
#          "H4K20me1" = 100, "H3K4me3" = 100, "H3K4me1" = 100,
#          "H3K36me3" = 100, "H3K4me2" = 100, "H3K9me3" = 100,
#          "H3K27me3" = 100)

lab <- c("expression" = 12, "H3K27ac" = 12, "H3K9ac" = 12,
         "H4K20me1" = 12, "H3K4me3" = 12, "H3K4me1" = 12,
         "H3K36me3" = 12, "H3K4me2" = 12, "H3K9me3" = 12,
         "H3K27me3" = 12)

lab2 <- names(lab)
kta1 <- ktab.data.frame(as.data.frame(w1), lab, tabnames = lab2)
statis1 <- statis(kta1, scannf = F)
plot(statis1)

df <- statis1$C.Co
df$lab <- rownames(df)
df <- cbind(df, do.call("rbind", strsplit(df$lab, split = "\\.")))

colnames(df)[5:6] <- c("tp", "mark")



ggplot(df, aes(C1, C2, size=tp)) +
  geom_point() +
  facet_wrap(~mark, ncol = 3)

###----------
### DISTATIS
###----------


# library(DistatisR)
# 
# expression <- as.matrix(dist(expression.matrix.scaled[, my.genes.safe]))
# H3K27ac <- as.matrix(dist(H3K27ac.matrix.scaled[, my.genes.safe]))
# H3K27me3 <- as.matrix(dist(H3K27me3.matrix.scaled[, my.genes.safe]))
# 
# mds<-cmdscale(dist(rbind(expression.matrix.scaled[, my.genes.safe],
#       H3K27me3.matrix.scaled[, my.genes.safe])))
# pca <- prcomp(rbind(expression.matrix.scaled[, my.genes.safe],
#                     H3K27me3.matrix.scaled[, my.genes.safe]))      
# 
# 
# A <- array(data = c(expression, H3K27ac, H3K27me3), c(12, 12, 3))
# d <- distatis(A)
# df <- d$res4Splus$F
# expression <- d$res4Splus$PartialF[,1:2,1]
# H3K27ac <- d$res4Splus$PartialF[,1:2,2]
# H3K27me3 <- d$res4Splus$PartialF[,1:2,3]
# 
# plot(expression, ylim = c(-.5,.5))
# points(H3K27ac, col = 2)
# points(H3K27me3, col = 3)



save.image(file="/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/figures/statis&distatis.RData")
