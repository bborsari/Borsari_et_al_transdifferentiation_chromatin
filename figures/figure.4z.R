.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


#************
# LIBRARIES *
#************

library(ggplot2)
library(reshape2)


palette <- c("no_peak" = "#000000",
             "mark_ant" = "#bdbdbd",
             "no_diff" = "#737373",
             "exp_ant" = "#403734",
             "stable" = "#fec44f",
             "not_correlated" = "#ec7014",
             "negatively_correlated" = "#993404")


hours <- c("1" = 0, "2" = 3, "3" = 6, "4" = 9,
           "5" = 12, "6" = 18, "7" = 24,
           "8" = 36, "9" = 48, "10" = 72,
           "11" = 120, "12" = 168)



#************
# FUNCTIONS *
#************

rescale <- function(x){
  return((x -min(x))/(max(x) - min(x)))
}


#********
# BEGIN *
#********

# 1. read expression matrix
expression.m <- read.table("H3K4me3/QN.merged/expression.matrix.tsv", h=T, sep="\t")

# 2. retrieve upregulated genes
metadata <- read.table("H3K4me3/QN.merged/metadata.tsv", h=T, sep="\t")
upreg.genes <- rownames(metadata[metadata$final_class == "upregulation", ])

# 3. subset expression matrix for upregulated genes
expression.m <- expression.m[rownames(expression.m) %in% upreg.genes, ]

# 4. keep upregulated genes in cluster 2
cluster2 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.2.txt", 
                       h=F, sep="\t")
expression.m <- expression.m[rownames(expression.m) %in% cluster2$V1, ]

# 5. order expression matrix from least to most expressed genes at 0h
expression.m <- expression.m[order(expression.m$H000, decreasing = F), ]

# 6. rescale expression profiles in range 0-1
expression.m.rescaled <- t(apply(expression.m, 1, rescale))

# 7. retrieve 
# 7.1. lowest expressed genes
lowest.expression.m.rescaled <- expression.m.rescaled[1:111, ]
lowest.expression.m <- expression.m[1:111, ]
lowest.expression.m$type2 <- "lowest"
lowest.expression.m.rescaled.melt <- melt(lowest.expression.m.rescaled)
lowest.expression.m.rescaled.melt$type <- "expression"
lowest.expression.m.rescaled.melt$type2 <- "lowest"

# 7.2. highest expressed genes
highest.expression.m.rescaled <- expression.m.rescaled[998:1108, ]
highest.expression.m <- expression.m[998:1108, ]
highest.expression.m$type2 <- "highest"
highest.expression.m.rescaled.melt <- melt(highest.expression.m.rescaled)
highest.expression.m.rescaled.melt$type <- "expression"
highest.expression.m.rescaled.melt$type2 <- "highest"

highest.lowest.expression.m <- rbind(highest.expression.m, lowest.expression.m)
highest.lowest.expression.m.melt <- melt(highest.lowest.expression.m)

highest.lowest.expression.m.melt$type2 <- factor(highest.lowest.expression.m.melt$type2,
                                                 levels = c("lowest", "highest"))

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.4z2.pdf",
    height = 3, width = 3)
ggplot(highest.lowest.expression.m.melt[highest.lowest.expression.m.melt$variable %in% c("H000"), ], 
       aes(x=type2, y=value)) +
  geom_jitter() +
  theme_bw() +
  theme(axis.title.y = element_text(size=13),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=15, vjust = .5),
        axis.text.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15, angle = 0),
        strip.background = element_blank(),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=15),
        legend.title = element_blank()) +
  ylab("expression at 0h - log2(TPM+1)")
dev.off()
  
  

# 8. source figure.4y.R
source("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/figures/figure.4y.R")



# 9. retrieve groups for lowest and highest expressed genes at 0h
# 9.1. lowest
df.barplot.lowest <- data.frame(stringsAsFactors = F)
for ( i in marks) {
  tmp.df <- as.data.frame(t(table(all.marks.groups[rownames(all.marks.groups) 
                                                   %in% rownames(lowest.expression.m) &
                                                     all.marks.groups[, i] != "positively_correlated", i])/
                              sum(table(all.marks.groups[rownames(all.marks.groups) 
                                                         %in% rownames(lowest.expression.m) &
                                                           all.marks.groups[, i] != "positively_correlated", i]))))
  tmp.df$mark <- i
  df.barplot.lowest <- rbind(df.barplot.lowest,
                             tmp.df)
}
df.barplot.lowest$Var1 <- "lowest"
df.barplot.lowest$mark <- factor(df.barplot.lowest$mark, levels = c("H3K4me1",
                                                                    "H3K4me2",
                                                                    "H3K27ac",
                                                                    "H3K9ac",
                                                                    "H3K4me3",
                                                                    "H3K36me3",
                                                                    "H4K20me1"))
df.barplot.lowest$Var2 <- factor(df.barplot.lowest$Var2, levels = c("negatively_correlated",
                                                                    "not_correlated",
                                                                    "stable",
                                                                    "mark_ant",
                                                                    "no_diff",
                                                                    "exp_ant",
                                                                    "no_peak"))

# 9.2. highest
df.barplot.highest <- data.frame(stringsAsFactors = F)
for ( i in marks) {
  tmp.df <- as.data.frame(t(table(all.marks.groups[rownames(all.marks.groups) 
                                                   %in% rownames(highest.expression.m) &
                                                     all.marks.groups[, i] != "positively_correlated", i])/
                              sum(table(all.marks.groups[rownames(all.marks.groups) 
                                                         %in% rownames(highest.expression.m) &
                                                           all.marks.groups[, i] != "positively_correlated", i]))))
  tmp.df$mark <- i
  df.barplot.highest <- rbind(df.barplot.highest,
                              tmp.df)
}
df.barplot.highest$Var1 <- "highest"
df.barplot.highest$mark <- factor(df.barplot.highest$mark, levels = c("H3K4me1",
                                                                      "H3K4me2",
                                                                      "H3K27ac",
                                                                      "H3K9ac",
                                                                      "H3K4me3",
                                                                      "H3K36me3",
                                                                      "H4K20me1"))
df.barplot.highest$Var2 <- factor(df.barplot.highest$Var2, levels = c("negatively_correlated",
                                                                      "not_correlated",
                                                                      "stable",
                                                                      "mark_ant",
                                                                      "no_diff",
                                                                      "exp_ant",
                                                                      "no_peak"))
df.barplot.highest <- df.barplot.highest[complete.cases(df.barplot.highest), ]

df.barplot <- rbind(df.barplot.lowest, df.barplot.highest)
df.barplot$Var1 <- factor(df.barplot$Var1, levels = c("lowest", "highest"))


# 10. plot
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.4z.pdf", 
    height = 6, width = 10)
ggplot(df.barplot, 
       aes(x=mark, y=Freq*100, fill=Var2)) +
  geom_bar(stat="identity", color = "white", alpha = .9) +
  scale_fill_manual(values = palette) +
  facet_grid(~Var1) +
  theme_bw() +
  theme(axis.title.y = element_text(size=11),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=15, angle = 30, vjust = .5),
        axis.text.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15, angle = 0),
        strip.background = element_blank(),
        plot.title = element_text(size = 22),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=15),
        legend.title = element_blank()) +
  ylab("% of genes")
dev.off()


write.table(as.data.frame(rownames(highest.expression.m)),
            file="~/public_html/Borsari_et_al_transdifferentiation_chromatin/GO.analysis/upregulated.2.highest.txt",
            quote=F, row.names = F, col.names = F, sep="\t")

write.table(as.data.frame(rownames(lowest.expression.m)),
            file="~/public_html/Borsari_et_al_transdifferentiation_chromatin/GO.analysis/upregulated.2.lowest.txt",
            quote=F, row.names = F, col.names = F, sep="\t")
