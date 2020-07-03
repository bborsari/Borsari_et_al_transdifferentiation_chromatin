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
  description = "\nPlot HMM emission matrix."
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options



#************
# LIBRARIES *
#************


library(pheatmap)
library(dplyr)
library(tidyr)
library(scales)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(viridis)
library(cowplot)
library(ggplotify)




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
                          rep("stably\ne.", length(stably.expressed.genes)),
                          rep("not\ne.", length(not.expressed.genes))))

# 7. set HMM wd
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/HMM/marks/")


# 8. read matrix of states' sequence for each gene
m <- read.table(paste0("HMM.", opt$n_states, ".gene.matrix.tsv"), h=T, sep="\t")


# 9. count, for each gene, how many states are found
m$states <- apply(m, 1, function(x){length(table(x))})
m$gene_id <- rownames(m)


# 10. add group info (not exp/stably exp/DE)
m <- merge(m, x, by = "gene_id")
m$group <- factor(m$group, levels = c("not\ne.", "stably\ne.", "DE"))
m$group2 <- "all"


# 11. summarize number of variable genes per group and in total
y <- melt(table(m[, c("states", "group")]))
z <- melt(table(m[, c("states", "group2")]))
colnames(z)[2] <- "group"
y <- rbind(y, z)


# 12. make plot
pdf(paste0("~/public_html/Borsari_et_al_transdifferentiation_chromatin/HMM/n_states/",
           opt$n_states, "/barplot.variable.genes.pdf"), width = 5, height = 4)
ggplot(y, aes(x=group, y=value, label = value, fill = as.factor(states))) +
  geom_bar(stat="identity", position = "fill") +
  labs(fill = "# of states\nin sequence") +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_brewer()
dev.off()