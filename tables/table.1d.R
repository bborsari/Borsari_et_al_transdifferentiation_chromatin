.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")



#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(reshape2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(rlist)
library(scales)
library(xtable)


#********
# BEGIN *
#********

# 0. source table.1g
source("/nfs/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/tables/table.1g.R")


# 1. set working directory
setwd("/nfs/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


# 2. read, for each mark, the ".6.groups.tsv" dataframe
marks <- c("H3K4me1", "H3K4me2", "H3K9ac", "H3K27ac", "H3K4me3",
           "H3K36me3", "H4K20me1", "H3K9me3", "H3K27me3")

lom <- list()

for ( i in 1:9 ) {
  
  tmp <- read.table(paste0(marks[i], "/QN.merged/", marks[i], ".6.groups.tsv"), 
                    h=T, sep="\t", 
                    stringsAsFactors = F)
  tmp$group <- gsub("peak_not_TSS", "no_peak", tmp$group) # merge "peak_not_TSS" and "no_peak"
  lom[[i]] <- tmp
  
}

rm(tmp)


# 3. build table of the 5 groups
m <- c() 

for (i in 1:9) {
  
  # 3.1. % of differentially marked genes divided into positively, not, and negatively c genes
  x <- table(lom[[i]]$group)[c("positively_correlated", "not_correlated", "negatively_correlated")]
  x <- paste0( x, " (", round( ((x/sum(x))*100), 2 ), "%)" )
  
  # 3.2. % of stably marked and not marked genes
  y <- table(lom[[i]]$group)[c("stable", "no_peak")]
  y <- paste0( y, " (", round( ((y/8030)*100), 2 ), "%)" )
  
  z <- as.character(c(x, y))
  names(z) <- c("positively_correlated", "not_correlated", "negatively_correlated",
                "no_peak", "stable")
  
  # 3.3. add x and y
  m <- rbind(m, z)
  
}

rownames(m) <- marks


# 4. merge m and m2
stopifnot(identical(rownames(m), m2$mark))
m <- cbind(m2, m[, 1:3])


# 5. save table
pdf("~bborsari/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/table.1d.pdf", 
    width=12, height = 4)
grid.table(m)
dev.off()


# 6. print table in latex format
print(xtable(m, caption = "DE genes (8030)"), include.rownames=F)

