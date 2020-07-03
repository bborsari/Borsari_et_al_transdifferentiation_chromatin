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



#************
# FUNCTIONS *
#************

rescale <- function(x){
  return((x -min(x))/(max(x) - min(x)))
}


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


# 8. define color palette
palette <- rev(rainbow(opt$n_states))
names(palette) <- z


# 9. make plot
pdf(paste0("~/public_html/Borsari_et_al_transdifferentiation_chromatin/HMM/n_states/",
           opt$n_states, "/emission.matrix.pdf"))
pheatmap(m[, 1:9], cluster_cols = F, cluster_rows = F,
         cellheight = 20, cellwidth = 20,
         border_color = NA,
         main = "emission matrix",
         color = c('#fff7f3','#fde0dd','#fcc5c0',
                   '#fa9fb5','#f768a1','#dd3497',
                   '#ae017e','#7a0177','#49006a'),
         labels_row = z,
         annotation_row = m[, c("states"), drop=F],
         annotation_colors = list(states = palette))
dev.off()
