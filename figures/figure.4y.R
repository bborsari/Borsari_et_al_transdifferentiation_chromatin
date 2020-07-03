.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(pheatmap)


#************
# FUNCTIONS *
#************

hours <- c("1" = 0, "2" = 3, "3" = 6, "4" = 9,
           "5" = 12, "6" = 18, "7" = 24,
           "8" = 36, "9" = 48, "10" = 72,
           "11" = 120, "12" = 168)


my.function1 <- function(f1.degree, f1.expression.mark.m, f1.integer.vector) {
  
  f1.my.list <- list(c(1,2,6), c(1,3,7), c(1,4,8), c(1,5,9))
  names(f1.my.list) <- c("perc_25", "perc_50", "perc_75", "perc_100")
  
  f1.expression.mark.m.subset <- f1.expression.mark.m[, f1.my.list[[f1.degree]]]
  colnames(f1.expression.mark.m.subset) <- c("gene_id", "expression", "mark")
  
  f1.integer.vector.subset <- f1.integer.vector[rownames(f1.integer.vector) 
                                                %in% f1.expression.mark.m.subset$gene_id, ]
  
  f1.exp.ant.df <- f1.expression.mark.m.subset[f1.expression.mark.m.subset$expression < f1.expression.mark.m.subset$mark, ]
  f1.exp.ant.shifts <- c()
  
  for ( i in 1:nrow(f1.exp.ant.df) ){
    
    f1.my.gene <- f1.exp.ant.df[i, "gene_id"]
    f1.my.exp.tp <- f1.exp.ant.df[i, 2]
    f1.my.alignment <- f1.integer.vector.subset[f1.my.gene, f1.my.exp.tp]
    f1.my.shift <- as.integer(hours[f1.my.exp.tp + f1.my.alignment]) - as.integer(hours[f1.my.exp.tp])
    f1.exp.ant.shifts <- c(f1.exp.ant.shifts, f1.my.shift)
    
  }
  
  f1.mark.ant.df <- f1.expression.mark.m.subset[f1.expression.mark.m.subset$expression > f1.expression.mark.m.subset$mark, ]
  f1.mark.ant.shifts <- c()
  for ( i in 1:nrow(f1.mark.ant.df) ){
    
    f1.my.gene <- f1.mark.ant.df[i, "gene_id"]
    f1.my.exp.tp <- f1.mark.ant.df[i, 2]
    f1.my.alignment <- f1.integer.vector.subset[f1.my.gene, f1.my.exp.tp]
    f1.my.shift <- as.integer(hours[f1.my.exp.tp + f1.my.alignment]) - as.integer(hours[f1.my.exp.tp])
    f1.mark.ant.shifts <- c(f1.mark.ant.shifts, f1.my.shift)
    
  }
  
  
  f1.no.diff.df <- f1.expression.mark.m.subset[f1.expression.mark.m.subset$expression == f1.expression.mark.m.subset$mark, ]
  f1.no.diff.shifts <- c()
  for ( i in 1:nrow(f1.no.diff.df) ){
    
    f1.my.gene <- f1.no.diff.df[i, "gene_id"]
    f1.my.exp.tp <- f1.no.diff.df[i, 2]
    f1.my.alignment <- f1.integer.vector.subset[f1.my.gene, f1.my.exp.tp]
    f1.my.shift <- as.integer(hours[f1.my.exp.tp + f1.my.alignment]) - as.integer(hours[f1.my.exp.tp])
    f1.no.diff.shifts <- c(f1.no.diff.shifts, f1.my.shift)
    
  }
  
  f1.my.out <- data.frame(shift = c(f1.exp.ant.shifts, f1.mark.ant.shifts, f1.no.diff.shifts),
                          group = c(rep("exp_ant", length(f1.exp.ant.shifts)),
                                    rep("mark_ant", length(f1.mark.ant.shifts)),
                                    rep("no_diff", length(f1.no.diff.shifts))),
                          degree = rep(f1.degree, nrow(f1.integer.vector.subset)),
                          gene_id = c(f1.exp.ant.df$gene_id, f1.mark.ant.df$gene_id, f1.no.diff.df$gene_id))
  
  return(f1.my.out)
  
}


my.function2 <- function(f2.mark) {
  
  # 1. read expression matrix
  f2.expression.m <- read.table(paste0(f2.mark, "/QN.merged/ant.del.analysis/increases/expression.25.50.75.100.tsv"),
                                h=T, sep="\t")
  f2.expression.m$gene_id <- rownames(f2.expression.m)
  
  # 2. read mark matrix
  f2.mark.m <- read.table(paste0(f2.mark, "/QN.merged/ant.del.analysis/increases/", f2.mark, ".25.50.75.100.tsv"),
                          h=T, sep="\t")
  f2.mark.m$gene_id <- rownames(f2.mark.m)
  
  # 3. merge mark and expression matrices
  f2.expression.mark.m <- merge(f2.expression.m, f2.mark.m, by="gene_id")
  
  # 4. keep only up-regulated and positively_correlated genes
  f2.upreg.pos.cc <- read.table(paste0(f2.mark, "/QN.merged/ant.del.analysis/increases/upregulation.positively_correlated.genes.txt"), 
                                h=F, sep="\t", stringsAsFactors = F)
  f2.upreg.pos.cc <- f2.upreg.pos.cc$V1
  f2.expression.mark.m <- f2.expression.mark.m[f2.expression.mark.m$gene_id %in% f2.upreg.pos.cc, ]
  
  # 5. read integer vector matrix
  f2.integer.vector <- read.table(paste0(f2.mark, "/QN.merged/ant.del.analysis/", f2.mark, ".integer.vector.tsv"), 
                                  h=T, sep = "\t")
  
  # 6. prepare df for plot
  f2.plot.df <- rbind(my.function1(f1.degree="perc_25", 
                                   f1.expression.mark.m = f2.expression.mark.m,
                                   f1.integer.vector = f2.integer.vector),
                      my.function1(f1.degree="perc_50", 
                                   f1.expression.mark.m = f2.expression.mark.m,
                                   f1.integer.vector = f2.integer.vector),
                      my.function1(f1.degree="perc_75", 
                                   f1.expression.mark.m = f2.expression.mark.m,
                                   f1.integer.vector = f2.integer.vector),
                      my.function1(f1.degree="perc_100", 
                                   f1.expression.mark.m = f2.expression.mark.m,
                                   f1.integer.vector = f2.integer.vector))
  
  f2.plot.df$mark <- f2.mark
  
  return(f2.plot.df)
  
}




#********
# BEGIN *
#********

# 1. read marks 6 groups dataframes
marks <- c("H3K4me1", "H3K4me2", "H3K4me3", "H3K27ac", "H3K9ac", "H3K36me3", "H4K20me1")

all.marks.groups.list <- list()

for ( i in 1:7 ) {
  
  tmp <- read.table(paste0(marks[i], "/QN.merged/", marks[i], ".6.groups.tsv"), h=T, sep="\t", stringsAsFactors = F)
  tmp$group <- gsub("peak_not_TSS", "no_peak", tmp$group)
  tmp$final_class <- gsub("regulation", "-\nregulated", tmp$final_class)
  all.marks.groups.list[[i]] <- tmp
  
}

rm(tmp)


# 2. check order of rownames
for (i in 2:7){
  
  stopifnot(identical(rownames(all.marks.groups.list[[1]]), 
                      rownames(all.marks.groups.list[[i]])))
}


# 3. prepare merged data.frame of groups across all marks
all.marks.groups <- data.frame(H3K4me1 = all.marks.groups.list[[1]]$group)

for (i in 2:7) {
  
  tmp <- data.frame(all.marks.groups.list[[i]]$group)
  colnames(tmp) <- marks[i]
  all.marks.groups <- cbind(all.marks.groups, tmp)
  
}

rm(tmp)
rownames(all.marks.groups) <- rownames(all.marks.groups.list[[1]])


# 3. keep upregulated genes only
metadata <- read.table("H3K4me3/QN.merged/metadata.tsv", h=T, sep="\t", stringsAsFactors = F)
upreg.genes <- metadata[metadata$final_class == "upregulation", ]
all.marks.groups <- all.marks.groups[rownames(all.marks.groups) %in% rownames(upreg.genes), ]


# 4. keep genes in cluster 2 
cluster2 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.2.txt")
all.marks.groups <- all.marks.groups[rownames(all.marks.groups) %in% cluster2$V1, ]


# 5. get group for 25% of upregulation 
my.plot.df <- data.frame(stringsAsFactors = F)

for ( i in marks ){
  
  my.plot.df <- rbind(my.plot.df, my.function2(f2.mark=i))
  
}

my.plot.df <- my.plot.df[my.plot.df$degree=="perc_25", ]


# 6. replace the category "positively_correlated" with exp_ant / mark_ant / no_diff

for ( i in marks ) {
  
  tmp.mark <- my.plot.df[my.plot.df$mark==i &
                           my.plot.df$gene_id %in% rownames(all.marks.groups), ]
  
  all.marks.groups[rownames(all.marks.groups) %in% tmp.mark$gene_id, "group"] <- tmp.mark$group
  
  all.marks.groups[, i] <- ifelse(is.na(all.marks.groups$group), 
                                     as.character(all.marks.groups[, i]), 
                                     as.character(all.marks.groups$group))
  all.marks.groups$group <- NULL
  
}





