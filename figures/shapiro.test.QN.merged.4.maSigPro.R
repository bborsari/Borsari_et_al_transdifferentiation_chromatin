.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


#************
# LIBRARIES *
#************

library(reshape2)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(xtable)


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

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


#************
# FUNCTIONS *
#************

f <- function(x) {
  if (diff(range(x)) == 0) NA else shapiro.test(x)$p.value
}





#********
# BEGIN *
#********

###------
### EXPRESSION
###-------
expression.R1 <- 
  read.table("expression/QN.merged/selected.genes.rep.2.after.QN.merged.tsv", 
             h=T, sep="\t")
expression.R1$gene_id <- rownames(expression.R1)

expression.R2 <- 
  read.table("expression/QN.merged/selected.genes.rep.3.after.QN.merged.tsv", 
             h=T, sep="\t")
expression.R2$gene_id <- rownames(expression.R2)

expression.R1.R2 <- merge(expression.R1, expression.R2, by = "gene_id")
rownames(expression.R1.R2) <- expression.R1.R2$gene_id
expression.R1.R2$gene_id <- NULL

pval.expression <- apply(expression.R1.R2, 1, f)
fdr.expression <- p.adjust(pval.expression, method="BH")

###----
### H3K4me3
###----
H3K4me3.R1 <- 
  read.table("H3K4me3/QN.merged/H3K4me3.R1.matrix.after.QN.merged.tsv", 
             h=T, sep="\t")
H3K4me3.R1$gene_id <- rownames(H3K4me3.R1)

H3K4me3.R2 <- 
  read.table("H3K4me3/QN.merged/H3K4me3.R2.matrix.after.QN.merged.tsv", 
             h=T, sep="\t")
H3K4me3.R2$gene_id <- rownames(H3K4me3.R2)

H3K4me3.R1.R2 <- merge(H3K4me3.R1, H3K4me3.R2, by = "gene_id")
H3K4me3.R1.R2 <- H3K4me3.R1.R2 %>% separate(gene_id, c("gene", "id"), "\\.")
rownames(H3K4me3.R1.R2) <- H3K4me3.R1.R2$gene
H3K4me3.R1.R2$gene <- NULL
H3K4me3.R1.R2$id <- NULL
H3K4me3.R1.R2 <- H3K4me3.R1.R2[rownames(expression.R1.R2), ]
stopifnot(identical(rownames(expression.R1.R2),
                    rownames(H3K4me3.R1.R2)))
colnames(H3K4me3.R1.R2) <- colnames(expression.R1.R2)

pval.H3K4me3 <- apply(H3K4me3.R1.R2, 1, f)
fdr.H3K4me3 <- p.adjust(pval.H3K4me3, method="BH")






###----
### H3K9ac
###----
H3K9ac.R1 <- 
  read.table("H3K9ac/QN.merged/H3K9ac.R1.matrix.after.QN.merged.tsv", 
             h=T, sep="\t")
H3K9ac.R1$gene_id <- rownames(H3K9ac.R1)

H3K9ac.R2 <- 
  read.table("H3K9ac/QN.merged/H3K9ac.R2.matrix.after.QN.merged.tsv", 
             h=T, sep="\t")
H3K9ac.R2$gene_id <- rownames(H3K9ac.R2)

H3K9ac.R1.R2 <- merge(H3K9ac.R1, H3K9ac.R2, by = "gene_id")
H3K9ac.R1.R2 <- H3K9ac.R1.R2 %>% separate(gene_id, c("gene", "id"), "\\.")
rownames(H3K9ac.R1.R2) <- H3K9ac.R1.R2$gene
H3K9ac.R1.R2$gene <- NULL
H3K9ac.R1.R2$id <- NULL
H3K9ac.R1.R2 <- H3K9ac.R1.R2[rownames(expression.R1.R2), ]
stopifnot(identical(rownames(expression.R1.R2),
                    rownames(H3K9ac.R1.R2)))
colnames(H3K9ac.R1.R2) <- colnames(expression.R1.R2)


pval.H3K9ac <- apply(H3K9ac.R1.R2, 1, f)
fdr.H3K9ac <- p.adjust(pval.H3K9ac, method="BH")



###----
### H3K27ac
###----
H3K27ac.R1 <- 
  read.table("H3K27ac/QN.merged/H3K27ac.R1.matrix.after.QN.merged.tsv", 
             h=T, sep="\t")
H3K27ac.R1$gene_id <- rownames(H3K27ac.R1)

H3K27ac.R2 <- 
  read.table("H3K27ac/QN.merged/H3K27ac.R2.matrix.after.QN.merged.tsv", 
             h=T, sep="\t")
H3K27ac.R2$gene_id <- rownames(H3K27ac.R2)

H3K27ac.R1.R2 <- merge(H3K27ac.R1, H3K27ac.R2, by = "gene_id")
H3K27ac.R1.R2 <- H3K27ac.R1.R2 %>% separate(gene_id, c("gene", "id"), "\\.")
rownames(H3K27ac.R1.R2) <- H3K27ac.R1.R2$gene
H3K27ac.R1.R2$gene <- NULL
H3K27ac.R1.R2$id <- NULL
H3K27ac.R1.R2 <- H3K27ac.R1.R2[rownames(expression.R1.R2), ]
stopifnot(identical(rownames(expression.R1.R2),
                    rownames(H3K27ac.R1.R2)))
colnames(H3K27ac.R1.R2) <- colnames(expression.R1.R2)


pval.H3K27ac <- apply(H3K27ac.R1.R2, 1, f)
fdr.H3K27ac <- p.adjust(pval.H3K27ac, method="BH")






###----
### H3K4me1
###----
H3K4me1.R1 <- 
  read.table("H3K4me1/QN.merged/H3K4me1.R1.matrix.after.QN.merged.tsv", 
             h=T, sep="\t")
H3K4me1.R1$gene_id <- rownames(H3K4me1.R1)

H3K4me1.R2 <- 
  read.table("H3K4me1/QN.merged/H3K4me1.R2.matrix.after.QN.merged.tsv", 
             h=T, sep="\t")
H3K4me1.R2$gene_id <- rownames(H3K4me1.R2)

H3K4me1.R1.R2 <- merge(H3K4me1.R1, H3K4me1.R2, by = "gene_id")
H3K4me1.R1.R2 <- H3K4me1.R1.R2 %>% separate(gene_id, c("gene", "id"), "\\.")
rownames(H3K4me1.R1.R2) <- H3K4me1.R1.R2$gene
H3K4me1.R1.R2$gene <- NULL
H3K4me1.R1.R2$id <- NULL
H3K4me1.R1.R2 <- H3K4me1.R1.R2[rownames(expression.R1.R2), ]
stopifnot(identical(rownames(expression.R1.R2),
                    rownames(H3K4me1.R1.R2)))
colnames(H3K4me1.R1.R2) <- colnames(expression.R1.R2)


pval.H3K4me1 <- apply(H3K4me1.R1.R2, 1, f)
fdr.H3K4me1 <- p.adjust(pval.H3K4me1, method="BH")



###----
### H3K4me2
###----
H3K4me2.R1 <- 
  read.table("H3K4me2/QN.merged/H3K4me2.R1.matrix.after.QN.merged.tsv", 
             h=T, sep="\t")
H3K4me2.R1$gene_id <- rownames(H3K4me2.R1)

H3K4me2.R2 <- 
  read.table("H3K4me2/QN.merged/H3K4me2.R2.matrix.after.QN.merged.tsv", 
             h=T, sep="\t")
H3K4me2.R2$gene_id <- rownames(H3K4me2.R2)

H3K4me2.R1.R2 <- merge(H3K4me2.R1, H3K4me2.R2, by = "gene_id")
H3K4me2.R1.R2 <- H3K4me2.R1.R2 %>% separate(gene_id, c("gene", "id"), "\\.")
rownames(H3K4me2.R1.R2) <- H3K4me2.R1.R2$gene
H3K4me2.R1.R2$gene <- NULL
H3K4me2.R1.R2$id <- NULL
H3K4me2.R1.R2 <- H3K4me2.R1.R2[rownames(expression.R1.R2), ]
stopifnot(identical(rownames(expression.R1.R2),
                    rownames(H3K4me2.R1.R2)))
colnames(H3K4me2.R1.R2) <- colnames(expression.R1.R2)


pval.H3K4me2 <- apply(H3K4me2.R1.R2, 1, f)
fdr.H3K4me2 <- p.adjust(pval.H3K4me2, method="BH")



###----
### H3K27me3
###----
H3K27me3.R1 <- 
  read.table("H3K27me3/QN.merged/H3K27me3.R1.matrix.after.QN.merged.tsv", 
             h=T, sep="\t")
H3K27me3.R1$gene_id <- rownames(H3K27me3.R1)

H3K27me3.R2 <- 
  read.table("H3K27me3/QN.merged/H3K27me3.R2.matrix.after.QN.merged.tsv", 
             h=T, sep="\t")
H3K27me3.R2$gene_id <- rownames(H3K27me3.R2)

H3K27me3.R1.R2 <- merge(H3K27me3.R1, H3K27me3.R2, by = "gene_id")
H3K27me3.R1.R2 <- H3K27me3.R1.R2 %>% separate(gene_id, c("gene", "id"), "\\.")
rownames(H3K27me3.R1.R2) <- H3K27me3.R1.R2$gene
H3K27me3.R1.R2$gene <- NULL
H3K27me3.R1.R2$id <- NULL
H3K27me3.R1.R2 <- H3K27me3.R1.R2[rownames(expression.R1.R2), ]
stopifnot(identical(rownames(expression.R1.R2),
                    rownames(H3K27me3.R1.R2)))
colnames(H3K27me3.R1.R2) <- colnames(expression.R1.R2)

pval.H3K27me3 <- apply(H3K27me3.R1.R2, 1, f)
fdr.H3K27me3 <- p.adjust(pval.H3K27me3, method="BH")



###----
### H3K36me3
###----
H3K36me3.R1 <- 
  read.table("H3K36me3/QN.merged/H3K36me3.R1.matrix.after.QN.merged.tsv", 
             h=T, sep="\t")
H3K36me3.R1$gene_id <- rownames(H3K36me3.R1)

H3K36me3.R2 <- 
  read.table("H3K36me3/QN.merged/H3K36me3.R2.matrix.after.QN.merged.tsv", 
             h=T, sep="\t")
H3K36me3.R2$gene_id <- rownames(H3K36me3.R2)

H3K36me3.R1.R2 <- merge(H3K36me3.R1, H3K36me3.R2, by = "gene_id")
H3K36me3.R1.R2 <- H3K36me3.R1.R2 %>% separate(gene_id, c("gene", "id"), "\\.")
rownames(H3K36me3.R1.R2) <- H3K36me3.R1.R2$gene
H3K36me3.R1.R2$gene <- NULL
H3K36me3.R1.R2$id <- NULL
H3K36me3.R1.R2 <- H3K36me3.R1.R2[rownames(expression.R1.R2), ]
stopifnot(identical(rownames(expression.R1.R2),
                    rownames(H3K36me3.R1.R2)))
colnames(H3K36me3.R1.R2) <- colnames(expression.R1.R2)


pval.H3K36me3 <- apply(H3K36me3.R1.R2, 1, f)
fdr.H3K36me3 <- p.adjust(pval.H3K36me3, method="BH")



###----
### H3K9me3
###----
H3K9me3.R1 <- 
  read.table("H3K9me3/QN.merged/H3K9me3.R1.matrix.after.QN.merged.tsv", 
             h=T, sep="\t")
H3K9me3.R1$gene_id <- rownames(H3K9me3.R1)

H3K9me3.R2 <- 
  read.table("H3K9me3/QN.merged/H3K9me3.R2.matrix.after.QN.merged.tsv", 
             h=T, sep="\t")
H3K9me3.R2$gene_id <- rownames(H3K9me3.R2)

H3K9me3.R1.R2 <- merge(H3K9me3.R1, H3K9me3.R2, by = "gene_id")
H3K9me3.R1.R2 <- H3K9me3.R1.R2 %>% separate(gene_id, c("gene", "id"), "\\.")
rownames(H3K9me3.R1.R2) <- H3K9me3.R1.R2$gene
H3K9me3.R1.R2$gene <- NULL
H3K9me3.R1.R2$id <- NULL
H3K9me3.R1.R2 <- H3K9me3.R1.R2[rownames(expression.R1.R2), ]
stopifnot(identical(rownames(expression.R1.R2),
                    rownames(H3K9me3.R1.R2)))
colnames(H3K9me3.R1.R2) <- colnames(expression.R1.R2)


pval.H3K9me3 <- apply(H3K9me3.R1.R2, 1, f)
fdr.H3K9me3 <- p.adjust(pval.H3K9me3, method="BH")



###----
### H4K20me1
###----
H4K20me1.R1 <- 
  read.table("H4K20me1/QN.merged/H4K20me1.R1.matrix.after.QN.merged.tsv", 
             h=T, sep="\t")
H4K20me1.R1$gene_id <- rownames(H4K20me1.R1)

H4K20me1.R2 <- 
  read.table("H4K20me1/QN.merged/H4K20me1.R2.matrix.after.QN.merged.tsv", 
             h=T, sep="\t")
H4K20me1.R2$gene_id <- rownames(H4K20me1.R2)

H4K20me1.R1.R2 <- merge(H4K20me1.R1, H4K20me1.R2, by = "gene_id")
H4K20me1.R1.R2 <- H4K20me1.R1.R2 %>% separate(gene_id, c("gene", "id"), "\\.")
rownames(H4K20me1.R1.R2) <- H4K20me1.R1.R2$gene
H4K20me1.R1.R2$gene <- NULL
H4K20me1.R1.R2$id <- NULL
H4K20me1.R1.R2 <- H4K20me1.R1.R2[rownames(expression.R1.R2), ]
stopifnot(identical(rownames(expression.R1.R2),
                    rownames(H4K20me1.R1.R2)))
colnames(H4K20me1.R1.R2) <- colnames(expression.R1.R2)


pval.H4K20me1 <- apply(H4K20me1.R1.R2, 1, f)
fdr.H4K20me1 <- p.adjust(pval.H4K20me1, method="BH")



###-------
### all marks
###-------


# distribution of FDR values
fdr.all.marks <- data.frame(fdr = c(fdr.expression,
                                    fdr.H3K4me3,
                                    fdr.H3K9ac,
                                    fdr.H3K27ac,
                                    fdr.H3K4me1,
                                    fdr.H3K4me2,
                                    fdr.H3K27me3,
                                    fdr.H3K36me3,
                                    fdr.H3K9me3,
                                    fdr.H4K20me1),
                            type = rep(c("expression", "H3K4me3", "H3K9ac",
                                         "H3K27ac", "H3K4me1",
                                         "H3K4me2", "H3K27me3",
                                         "H3K36me3", "H3K9me3", "H4K20me1"), each=12248))


fdr.all.marks$type <- factor(fdr.all.marks$type, levels = c("expression", "H3K27ac", "H3K9ac", "H3K4me3", 
                                                            "H3K4me1", "H3K4me2", "H4K20me1", 
                                                            "H3K36me3", "H3K9me3", "H3K27me3"))


for (el in c("expression", "H3K4me3", "H3K9ac",
             "H3K27ac", "H3K4me1",
             "H3K4me2", "H3K27me3",
             "H3K36me3", "H3K9me3", "H4K20me1")) {
  
  print(el)
  print (sum(is.na(fdr.all.marks[fdr.all.marks$type==el,])))
  
}


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/QN.merged.FDR.distribution.pdf", width=10, height=4.5)
ggplot(fdr.all.marks, aes(y=-log10(fdr), x=type, fill=type)) + 
  geom_hline(yintercept = 2, color="red", linetype="dashed") +
  geom_violin(alpha=.5, width=.9) +
  geom_boxplot(width=0.2, alpha = .6, outlier.shape = NA) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15, angle=30, vjust=.5),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ylab("-log10(FDR) - shapiro-wilk test") +
  scale_fill_manual(values = palette) +
  guides(fill=F)
dev.off()


# table of n. of genes passing shapiro-wilk test
fdr.table.all.marks <- rbind(table(fdr.expression < 0.01),
                             table(fdr.H3K4me3 < 0.01),
                             table(fdr.H3K9ac < 0.01),
                             table(fdr.H3K27ac < 0.01),
                             table(fdr.H3K4me1 < 0.01),
                             table(fdr.H3K4me2 < 0.01),
                             table(fdr.H3K27me3 < 0.01),
                             table(fdr.H3K36me3 < 0.01),
                             table(fdr.H3K9me3 < 0.01),
                             table(fdr.H4K20me1 < 0.01))

rownames(fdr.table.all.marks) <- c("expression", "H3K4me3", "H3K9ac",
                                   "H3K27ac", "H3K4me1",
                                   "H3K4me2", "H3K27me3", "H3K36me3",
                                   "H3K9me3", "H4K20me1")                        
colnames(fdr.table.all.marks) <- c("normal profile", "not normal profile")


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/QN.merged.FDR.table.pdf")
grid.table(fdr.table.all.marks)
# xtable(fdr.table.all.marks)
dev.off()