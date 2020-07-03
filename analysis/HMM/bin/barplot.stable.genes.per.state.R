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


# 8. set HMM wd
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/HMM/marks/")


# 9. read emission matrix
m <- read.table(paste0("HMM.", opt$n_states, ".response.matrix.tsv"), 
                h=T, sep="\t")

# 10. discard sd values
m <- m[, paste0("Re", 1:9, "..Intercept.")]


# 11. update colnames
colnames(m) <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3",
                 "H3K36me3", "H4K20me1", "H3K9me3", "H3K27me3")


# 12. sort states according to mean row value
m$mean <- apply(m, 1, mean)
m <- m[order(m$mean), ]
m$mean <- NULL


# 13. assign to each state a letter
z <- paste0(letters[1:opt$n_states], ":", gsub("St", "", rownames(m)))
names(z) <- rownames(m)


# 14. rescale marks values to range 0-1
m <- as.data.frame(apply(m, 2, rescale))
m$states <- z


# 15. read matrix of states' sequence for each gene
m2 <- read.table(paste0("HMM.", opt$n_states, ".gene.matrix.tsv"), 
                 h=T, sep="\t")


# 16. count, for each gene, how many and which states are found
m2$n_states <- apply(m2, 1, function(x){length(table(x))})
m2$states <- apply(m2[, 1:12], 1, function(x){paste(as.vector(unique(x)), collapse =";")})
m2$states2 <- z[paste0("St", m2$states)]
m2$gene_id <- rownames(m2)


# 17. add group info (not exp/stably exp/DE)
m2 <- merge(m2, x, by = "gene_id")
m2$group <- factor(m2$group, levels = c("not e.", "stably e.", "DE"))
m2$group2 <- "all"



# 18. define color palette
palette <- rev(rainbow(opt$n_states))
names(palette) <- z



# 19. make plots
lop <- list()


lop[[1]] <- ggplot(m2[m2$n_states < 2, ], aes(x=states2, fill=states2)) +
  geom_bar() +
  scale_fill_manual(values = palette) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.title = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = .5)) +
  labs(title = "all genes")


lop[[2]] <- ggplot(m2[m2$n_states < 2, ], aes(x=states2, fill=states2)) +
  geom_bar() +
  facet_grid(~group) +
  scale_fill_manual(values = palette) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 14),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.title = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = .5))
  



pdf(paste0("~/public_html/Borsari_et_al_transdifferentiation_chromatin/HMM/n_states/",
           opt$n_states, "/barplot.stable.genes.per.state.pdf"), width = 5, height = 10)
plot_grid(plotlist = lop, nrow = 4, ncol=1, scale = 0.9)
dev.off()



