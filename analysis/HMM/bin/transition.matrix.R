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

# 1. set HMM wd
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/HMM/marks/")


# 2. read emission matrix
m <- read.table(paste0("HMM.", opt$n_states, ".response.matrix.tsv"), 
                h=T, sep="\t")

# 3. discard sd values
m <- m[, paste0("Re", 1:9, "..Intercept.")]


# 4. update colnames
colnames(m) <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3",
                 "H3K36me3", "H4K20me1", "H3K9me3", "H3K27me3")


# 5. sort states according to mean row value
m$mean <- apply(m, 1, mean)
m <- m[order(m$mean), ]
m$mean <- NULL


# 6. assign to each state a letter
z <- paste0(letters[1:opt$n_states], ":", gsub("St", "", rownames(m)))
names(z) <- rownames(m)


# 7. rescale marks values to range 0-1
m <- as.data.frame(apply(m, 2, rescale))
m$states <- z


# 8. read transition matrix
m2 <- read.table(paste0("HMM.", opt$n_states, ".transition.matrix.tsv"), h=T, sep="\t")


# 9. update row- and colnames
colnames(m2) <- paste0("St", seq(1:opt$n_states))
rownames(m2) <- colnames(m2)


# 10. reorder states according to emission matrix
m2 <- m2[rownames(m), rownames(m)]
colnames(m2) <- z
rownames(m2) <- z


# 11. make plot
pdf(paste0("~/public_html/Borsari_et_al_transdifferentiation_chromatin/HMM/n_states/",
           opt$n_states, "/transition.matrix.pdf"))
pheatmap(m2, cluster_rows = F, cluster_cols = F,
         cellheight = 40, cellwidth = 40,
         display_numbers = T,
         number_color = "white",
         border_color = NA,
         main = "transition matrix")
dev.off()

