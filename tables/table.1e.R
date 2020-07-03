.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


#************
# LIBRARIES *
#************

library(gridExtra)
library(xtable)



#********
# BEGIN *
#********


# 1. set working directory
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


# 2. read table with stats
m <- read.table("silent.genes.no_peak.variable.stable.tsv", h=F, sep="\t")
colnames(m) <- c("mark", "no_peak", "variable", "stable")


# 3. reorder columns of m
m <- m[, c("mark", "no_peak", "stable", "variable")]


# 4. compute frequencies
m$no_peak_frac <- round((m$no_peak / 1552)*100, 2)
m$variable_frac <- round((m$variable / 1552)*100, 2)
m$stable_frac <- round((m$stable / 1552)*100, 2)


# 5. sanity check
apply(m[, 5:7], 1, sum)

m$no_peak <- paste0(m$no_peak, " (", m$no_peak_frac, "%)")
m$variable <- paste0(m$variable, " (", m$variable_frac, "%)")
m$stable <- paste0(m$stable, " (", m$stable_frac, "%)")


# 6. reorder marks
m$mark <- factor(m$mark, levels = c("H3K4me1", "H3K4me2", "H3K9ac", "H3K27ac", "H3K4me3",
                                    "H3K36me3", "H4K20me1", "H3K9me3", "H3K27me3"))
m <- m[order(m$mark), ]



# 7. save table
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/table.1e.pdf", 
    width = 10)
grid.table(m[, 1:4], rows = NULL)
dev.off()


# 8. table in latex format
print(xtable(m[, 1:4], caption = "not expressed genes - 1552"), include.rownames=FALSE)
