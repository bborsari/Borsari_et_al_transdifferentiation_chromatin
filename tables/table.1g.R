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


#********
# BEGIN *
#********


# 1. set working directory
setwd("/nfs/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


# 2. the marks we're analyzing
marks <- c("H3K4me1", "H3K4me2", "H3K9ac", "H3K27ac", "H3K4me3",
           "H3K36me3", "H4K20me1", "H3K9me3", "H3K27me3")

# 3. store, for each mark, the "6.groups.tsv" df
x <- list()

for ( i in 1:9 ) {
  
  tmp <- read.table(paste0(marks[i], "/QN.merged/", marks[i], ".6.groups.tsv"), 
                    h=T, sep="\t", stringsAsFactors = F)
  
  tmp$group <- gsub("peak_not_TSS", "no_peak", tmp$group) # merge peak_not_TSS and no_peak
  tmp$group <- gsub("not_correlated", "variable", tmp$group) # not_correlated is DM
  tmp$group <- gsub("negatively_correlated", "variable", tmp$group) # negatively_correlated is DM
  tmp$group <- gsub("positively_correlated", "variable", tmp$group) # positively_correlated is DM
  
  x[[i]] <- tmp
  
}

rm(tmp)


# 4. build table of frequency of not marked, stable and DM groups
m <- data.frame(mark = character(9),
                no_peak = integer(9),
                stable = integer(9),
                variable = integer(9),
                stringsAsFactors = F) 

for (i in 1:9) {
  
  # 4.1. compute frequency of the 3 groups
  y <- paste0(table(x[[i]]$group), " (", round(table(x[[i]]$group)/8030*100, 2), "%)" )
  
  # 4.2. add mark info
  y <- c( marks[i], y )
  names(y)[1] <- "mark"
  
  # 4.3. store vector
  m[i, ] <- y
  
}


# 5. reorder columns of m
m2 <- m[, c("mark", "no_peak", "stable", "variable")]


# 6. print table in latex format
print(xtable(m), include.rownames=F)
