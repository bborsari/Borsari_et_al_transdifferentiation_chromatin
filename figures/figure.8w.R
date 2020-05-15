.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")



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



#************
# FUNCTIONS *
#************

function1 <- function(mark) {
  
  df <- read.table(paste0(mark, "/QN.merged/ant.del.analysis/increases/", 
                          "expression.", mark, ".dynamics.tsv"),
                   h=T, sep="\t")
  
  df <- df[, c("group_25", "group_50", "group_75", "group_100")]
  
  df$gene_id <- rownames(df)
  
  df.melt <- melt(df, id.vars = "gene_id")
  df.melt$value <- factor(df.melt$value, levels = c("delayed", "concomitant", "anticipated"))
  
  df.melt$mark <- mark
  
  return(df.melt)
  
}


function2 <- function(df, m1, m2) {
  
  df.m1 <- df[df$mark == m1, ]
  df.m2 <- df[df$mark == m2, ]
  
  # create a merged dataframe on the common genes 
  # between mark1 and mark2
  df.m1.m2 <- merge(df.m1[, c("gene_id", "value")],
                    df.m2[, c("gene_id", "value")],
                    by = "gene_id")
  
  # update colnames
  colnames(df.m1.m2)[2:3] <- c(m1, m2)
  
  return(df.m1.m2)
  
}


function3 <- function(m, degree) {
  
  
  # 1. focus on groups defined considering the 
  # specified degree of up-regulation
  k <- m[m$variable==degree, ]
  
  # 2. compute overlap between all possible 
  # pairs of marks and groups (anticipated, concomitant, delayed)
  all.marks.matrix.overlap <- data.frame(stringsAsFactors = F)
  
  
  groups1 <- c("anticipated", "concomitant", "delayed")
  groups2 <- groups1
  
  for ( m1 in marks ) {
    
    overlap.vector <- c()
    
    for ( G1 in groups1 ) {
      
      for ( G2 in groups2 ) {
        
        for ( m2 in marks ) {
          
          print(paste(m1, G1, m2, G2)) # the combination of marks and groups currently analyzed
          
          # 3. obtain df of common genes between m1 and m2
          j <- function2(df = k, m1 = m1, m2 = m2)
          
          # 4. number of genes that are classified G1 for m1
          m1G1 <- nrow(j[j[, m1] == G1, ])
          
          # 5. number of genes classified G1 for m1, 
          # and classified G2 for m2
          m1G1.m2G2 <- nrow(j[j[, m2] == G2 & j[, m1] == G1, ])
          
          # 6. compute fraction of m1G1 over m1G1.m2G2
          overlap = (m1G1.m2G2 / m1G1)
          
          # 7. store overlap for the combination m1-G1, m2-G2
          overlap.vector <- c(overlap.vector, overlap)
          
        }
        
      }
      
    }
    
    
    # 8. compute matrix for m1 across all groups and marks
    matrix.overlap <- matrix(overlap.vector, ncol=21, nrow=3, byrow = T)
    
    
    # 9. add rownames and colnames for matrix.overlap
    rownames(matrix.overlap) <- paste(groups1, m1, sep="_")
    colnames(matrix.overlap) <- paste(rep(groups2, each = 7), rep(marks, 3), sep="_")
    
    
    # 10. store results for m1
    all.marks.matrix.overlap <- rbind(all.marks.matrix.overlap, matrix.overlap)
    
  }
  
  
  # 11. keep same order for rows and columns for matrix of all marks 
  all.marks.matrix.overlap.reordered  <- all.marks.matrix.overlap[colnames(all.marks.matrix.overlap), ]
  
  
  # 12. prepare annotation df for heatmap
  annotation.df <- data.frame(group = rep(groups1, each=7),
                              mark = rep(paste(marks, "   "), 3))
  rownames(annotation.df) <- rownames(all.marks.matrix.overlap.reordered)
  
  
  # 13. save the current matrix
  write.table(all.marks.matrix.overlap.reordered,
              paste0("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8w.",
                     degree,
                     ".tsv"),
              row.names = T, col.names = T, sep = "\t", quote=F)
  
  
  # 14. make heatmap
  p <- pheatmap(as.matrix(all.marks.matrix.overlap.reordered*100)[c(1:7, 15:21), ],
                color=c('white', '#fef5f4', '#fef2f1', "#fde0dd", '#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a'),
                breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
                cluster_rows = F,
                cluster_cols = F,
                border_color = NA,
                cellwidth = 15,
                cellheight = 15,
                annotation_row = annotation.df,
                annotation_col = annotation.df,
                annotation_colors = list(group = palette, mark = palette2),
                labels_row = annotation.df[c(1:7, 15:21), "mark"],
                show_colnames = F,
                annotation_legend = F,
                annotation_names_row = F,
                annotation_names_col = F,
                fontsize = 15,
                main = titles[degree])
  
  return(p)
  
  
}





#********
# BEGIN *
#********

# 1. set working directory
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


# 2. define palette for groups and marks
palette <- c("anticipated" = "#bdbdbd",
             "concomitant" = "#737373",
             "delayed" = "#403734")

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



# 3. the marks we're analyzing
marks <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3", "H3K36me3", "H4K20me1")


# 4. read dfs of all marks using function1
x <- data.frame(stringsAsFactors = F)
for ( i in 1:7) {
  
  x <- rbind(x, function1(mark = marks[i]))
  
}


# 5. titles for the 4 degrees of upregulation
titles <- c("group_25" = "25%",
            "group_50" = "50%",
            "group_75" = "75%",
            "group_100" = "100%")



# 6. make plots
# 6.1. 25% upregulation
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8w.25.pdf")
print(function3(m = x, degree = "group_25"))
dev.off()


# 6.2. 50% upregulation
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8w.50.pdf")
print(function3(m = x, degree = "group_50"))
dev.off()


# 6.3. 75% upregulation
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8w.75.pdf")
print(function3(m = x, degree = "group_75"))
dev.off()


# 6.4. 100% upregulation
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8w.100.pdf")
print(function3(m = x, degree = "group_100"))
dev.off()

