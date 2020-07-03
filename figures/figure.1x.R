.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")



#************
# LIBRARIES *
#************


library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(circlize)
library(seriation)
library(reshape2)
library(ggpubr)



#********
# BEGIN *
#********

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/")

# 1. read expression matrix
expression.matrix <- read.table("expression/QN.merged/selected.genes.rep.2.3.after.QN.merged.tsv",
                                h=T, sep="\t")

# 2. compute min, max and mean expression
expression.matrix$min <- apply(expression.matrix[, 1:12], 1, min)
expression.matrix$max <- apply(expression.matrix[, 1:12], 1, max)
expression.matrix$mean <- apply(expression.matrix[, 1:12], 1, mean)


# 3. retrieve flat genes for expression
flat.genes <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/genes/flat.genes.expression/flat.genes.txt",
                         h=F, sep="\t", stringsAsFactors = F)
flat.genes <- flat.genes$V1

# 4. retrieve flat genes for expression variable at least for one mark
flat.genes.variable.for.one.mark <- 
  read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/genes/flat.genes.expression/flat.genes.variable.for.one.mark.peaks.txt",
             h=F, sep="\t", stringsAsFactors = F)
flat.genes.variable.for.one.mark <- flat.genes.variable.for.one.mark$V1

# 5. retrieve flat genes for expression that are also flat for all marks
flat.genes.not.variable <- setdiff(flat.genes, flat.genes.variable.for.one.mark)

# 6. retrieve expression matrix for the two groups of genes
flat.genes.not.variable.e <- expression.matrix[rownames(expression.matrix) %in% flat.genes.not.variable, 
                                               c("min", "max", "mean")]
flat.genes.not.variable.e$group <- "exp_and_mark"

flat.genes.variable.for.one.mark.e <- expression.matrix[rownames(expression.matrix) %in% flat.genes.variable.for.one.mark, 
                                                        c("min", "max", "mean")]
flat.genes.variable.for.one.mark.e$group <- "exp_only"


# 7. genes for plot
m <- rbind(flat.genes.variable.for.one.mark.e, 
           flat.genes.not.variable.e)
m <- melt(m)

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1x.pdf",
    height = 4)
ggplot(m, aes(x=group, y=value)) +
  geom_boxplot() +
  facet_grid(~variable) +
  stat_compare_means() +
  theme_bw() +
  labs(title = "flat genes") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size=15, hjust = .5),
        axis.line = element_line(colour = "black")) +
  ylab("log2(TPM+1)") 
dev.off()
