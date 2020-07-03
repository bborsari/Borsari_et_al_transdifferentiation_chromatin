.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


#************
# LIBRARIES *
#************

library(gridExtra)

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")



#********
# BEGIN *
#********

m <- read.table("percentage.significant.genes.only.mark.silent.tsv", h=T, sep="\t")
m$percentage <- round(m$percentage*100, 2)
m <- m[order(m$n_genes, decreasing = T), ]


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/table.1c.pdf", 
    width = 10)
grid.table(m, rows = NULL)
dev.off()

print(xtable(m, caption = "silent genes (n = 1552)"), include.rownames=FALSE)
