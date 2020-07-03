.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(reshape2)
library(ggalluvial)
library(viridis)
library(ggridges)
library(pheatmap)
library(grid)
library(gtable)



#*************************
# PALETTES & TIME-POINTS *
#*************************

palette <- c("mark_ant" = "#bdbdbd",
             "no_diff" = "#737373",
             "exp_ant" = "#000000")

palette2 <- c("H3K9ac" = "#e7298a",
              "H3K27ac" = "#8e0152",
              "H3K4me3" = "#f1b6da",
              "H3K27me3" = "#253494",
              "H3K9me3" = "#41b6c4",
              "H3K36me3" = "#7fbc41",
              "H4K20me1" = "#276419",
              "H3K4me1" = "#ffbf00",
              "H3K4me2" = "#a67c00")
names(palette2) <- paste(names(palette2), "   ")

hours <- c("1" = 0, "2" = 3, "3" = 6, "4" = 9,
           "5" = 12, "6" = 18, "7" = 24,
           "8" = 36, "9" = 48, "10" = 72,
           "11" = 120, "12" = 168)



#************
# FUNCTIONS *
#************


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



my.function3 <- function(f3.mark1, f3.mark2, f3.plot.df) {
  
  f3.my.plot.df.mark1 <- f3.plot.df[f3.plot.df$mark==f3.mark1, ]
  f3.my.plot.df.mark2 <- f3.plot.df[f3.plot.df$mark==f3.mark2, ]
  
  # common genes
  f3.mark1.mark2.genes <- intersect(f3.my.plot.df.mark1$gene_id, f3.my.plot.df.mark2$gene_id)
  
  # create a merged dataframe on the common genes 
  f3.mark1.mark2.df <- merge(f3.my.plot.df.mark1[f3.my.plot.df.mark1$gene_id %in% f3.mark1.mark2.genes, 
                                                 c("group", "gene_id")],
                             f3.my.plot.df.mark2[f3.my.plot.df.mark2$gene_id %in% f3.mark1.mark2.genes, 
                                                 c("group", "gene_id")], by = "gene_id")
  # update colnames
  colnames(f3.mark1.mark2.df)[2:3] <- c(f3.mark1, f3.mark2)
  
  # find weird connections
  f3.my.linetype <- c()
  
  for ( i in 1:nrow(f3.mark1.mark2.df) ){
    
    if ( ( (f3.mark1.mark2.df[i, f3.mark1] == "no_diff") & 
           (f3.mark1.mark2.df[i, f3.mark2] == "mark_ant") ) |
         
         ( (f3.mark1.mark2.df[i, f3.mark1] == "exp_ant") & 
           (f3.mark1.mark2.df[i, f3.mark2] == "mark_ant") ) |
         
         ( (f3.mark1.mark2.df[i, f3.mark1] == "exp_ant") & 
           (f3.mark1.mark2.df[i, f3.mark2] == "no_diff") ) ) {
      
      f3.my.linetype <- c(f3.my.linetype, "weird")
      
    } else {
      
      f3.my.linetype <- c(f3.my.linetype, "OK")
      
    }
    
  }
  
  f3.mark1.mark2.df$my.linetype <- f3.my.linetype
  return(f3.mark1.mark2.df)
  
}



my.function4 <- function(f4.plot.df, f4.degree) {
  
  # 2. focus on groups defined considering 25% of up-regulation
  f4.test <- f4.plot.df[f4.plot.df$degree==f4.degree, ]
  
  # 3. compute overlap between all possible pairs
  f4.all.marks.matrix.overlap <- data.frame(stringsAsFactors = F)
  
  
  f4.groups1 <- c("exp_ant", "mark_ant", "no_diff")
  f4.groups2 <- f4.groups1
  f4.groups3 <- f4.groups2
  
  for ( f4.m1 in marks ) {
    
    f4.overlap.vector <- c()
    
    for ( f4.G1 in f4.groups1 ) {
      
      for ( f4.G2 in f4.groups2 ) {
        
        for ( f4.m2 in marks ) {
          
          print(paste(f4.m1, f4.G1, f4.m2, f4.G2))
          
          f4.mark1.mark2.df <- my.function3(f3.mark1 = f4.m1, f3.mark2 = f4.m2, f3.plot.df = f4.test)
          
          f4.mark1.mark2.genes <- f4.mark1.mark2.df[f4.mark1.mark2.df[, f4.m2] == f4.G2 & f4.mark1.mark2.df[, f4.m1] == f4.G1, "gene_id"]
          
          f4.expression.m <- read.table("H3K4me3/QN.merged/expression.matrix.tsv", h=T, sep="\t")
          f4.expression.m <- f4.expression.m[rownames(f4.expression.m) %in% f4.mark1.mark2.genes, ]
          f4.overlap <- median(f4.expression.m$H000)
          
          f4.overlap.vector <- c(f4.overlap.vector, f4.overlap)
          
          
        }
        
      }
      
    }
    
    f4.matrix.overlap <- matrix(f4.overlap.vector, ncol=21, nrow=3, byrow = T)
    
    f4.tmp.colnames <- c()
    for ( f4.G3 in f4.groups3 ) {
      
      f4.tmp.colnames <- c(f4.tmp.colnames, rep(f4.G3, 7))
      
    }
    
    colnames(f4.matrix.overlap) <- paste(f4.tmp.colnames, marks, sep="_")
    rownames(f4.matrix.overlap) <- paste(f4.groups1, f4.m1, sep="_")
    
    f4.all.marks.matrix.overlap <- rbind(f4.all.marks.matrix.overlap, f4.matrix.overlap)
    
  }
  
  f4.my.annotation.df <- data.frame(group = rep(c("exp_ant", "mark_ant", "no_diff"), 7),
                                    mark = rep(paste(marks, "   "), each=3))
  
  rownames(f4.my.annotation.df) <- paste(group = rep(c("exp_ant", "mark_ant", "no_diff"), 7),
                                          mark = rep(paste(marks), each=3), sep="_")
  
  f4.rows <- paste(rep(c("mark_ant", "no_diff", "exp_ant"), each=7),
                   rep(marks, 3), sep="_")
  
  
  # f4.rows <- names(sort(apply(f4.all.marks.matrix.overlap, 1, mean, na.rm=T)))
  # f4.cols <- names(sort(apply(f4.all.marks.matrix.overlap, 2, mean, na.rm=T)))
  
  f4.all.marks.matrix.overlap <- f4.all.marks.matrix.overlap[f4.rows, f4.rows]

  f4.p <- pheatmap(f4.all.marks.matrix.overlap[c(1:7, 15:21), ],
                   # color=c('#fff7f3', '#fde0dd', '#fcc5c0', '#fa9fb5', '#f768a1', '#dd3497', "#ae017e"),
                   # breaks = c(0, 1.5, 2.5, 3, 3.5, 4, 4.5, 5),
                   cluster_rows = F,
                   cluster_cols = F,
                   border_color = NA,
                   cellwidth = 15,
                   cellheight = 15,
                   annotation_row = f4.my.annotation.df,
                   # annotation_col = f4.my.annotation.df[f4.cols, c("group"), drop=F],
                   annotation_col = f4.my.annotation.df,
                   annotation_colors = list(group = palette,
                                            mark = palette2),
                   labels_row = f4.my.annotation.df[f4.rows, "mark"][c(1:7, 15:21)],
                   # labels_col = f4.my.annotation.df[f4.cols, "mark"],
                   # labels_col = f4.my.annotation.df[f4.rows, "mark"],
                   show_colnames = F,
                   annotation_legend = F,
                   annotation_names_row = F,
                   annotation_names_col = F,
                   fontsize = 15,
                   na_col = "#d9d9d9")
  
  return(f4.p)
  

  

}


#********
# BEGIN *
#********

marks <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", 
           "H3K4me3", "H3K36me3", "H4K20me1")

# 1. retrieve groups
my.plot.df <- data.frame(stringsAsFactors = F)

for ( i in marks ){
  
  my.plot.df <- rbind(my.plot.df, my.function2(f2.mark=i))
  
}


# 2. 25% upregulation
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.4w.25.pdf")
print(my.function4(f4.plot.df = my.plot.df,
                   f4.degree = "perc_25"))
dev.off()


# # 3. 50% upregulation
# pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.4w.50.pdf")
# print(my.function4(f4.plot.df = my.plot.df,
#                    f4.degree = "perc_50"))
# dev.off()
# 
# 
# # 4. 75% upregulation
# pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.4w.75.pdf")
# print(my.function4(f4.plot.df = my.plot.df,
#                    f4.degree = "perc_75"))
# dev.off()
# 
# 
# # 5. 100% upregulation
# pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.4w.100.pdf")
# print(my.function4(f4.plot.df = my.plot.df,
#                    f4.degree = "perc_100"))
# dev.off()
