.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")

# palette <- c("positively_correlated" = "#DEA450",
#              "no_peak" = "#9B461F",
#              "negatively_correlated" = "#a095a0",
#              "stable" = "#5d6f62",
#              "not_correlated" = "#42354C")

palette <- c("positively_correlated" = "#5B9F80",
              "no_peak" = "#403734",
              "negatively_correlated" = "#993404",
              "stable" = "#fec44f",
              "not_correlated" = "#ec7014")



#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(pheatmap)


#********
# BEGIN *
#********

# 1. read marks 6 groups dataframes
marks <- c("H3K27ac", "H3K9ac", "H4K20me1", "H3K4me3", "H3K4me1",
           "H3K36me3", "H3K4me2", "H3K9me3", "H3K27me3")

all.marks.groups.list <- list()

for ( i in 1:9 ) {
  
  tmp <- read.table(paste0(marks[i], "/QN.merged/", marks[i], ".6.groups.tsv"), h=T, sep="\t", stringsAsFactors = F)
  tmp$group <- gsub("peak_not_TSS", "no_peak", tmp$group)
  # tmp$group <- gsub("no_peak", "unmarked", tmp$group)
  
  tmp$final_class <- gsub("regulation", "-\nregulated", tmp$final_class)
  all.marks.groups.list[[i]] <- tmp
  
}

rm(tmp)


# 2. check order of rownames
for (i in 2:9){
  
  stopifnot(identical(rownames(all.marks.groups.list[[1]]), 
                      rownames(all.marks.groups.list[[i]])))
}


# 3. prepare merged data.frame of groups across all marks
all.marks.groups <- data.frame(H3K27ac = all.marks.groups.list[[1]]$group)

for (i in 2:9) {
  
  tmp <- data.frame(all.marks.groups.list[[i]]$group)
  colnames(tmp) <- marks[i]
  all.marks.groups <- cbind(all.marks.groups, tmp)
  
}

rm(tmp)
rownames(all.marks.groups) <- rownames(all.marks.groups.list[[1]])



# 4. perform hypergeometric test

groups1 <- c("negatively_correlated", "no_peak", "not_correlated", "positively_correlated", "stable")
groups2 <- groups1
groups3 <- groups2

all.marks.matrix.pvalues <- data.frame(stringsAsFactors = F)
all.marks.matrix.fdr <- data.frame(stringsAsFactors = F)

for ( m1 in marks ) {
  
  pvalues <- c()

  for ( G1 in groups1 ) {

    for ( G2 in groups2 ) {

      for ( m2 in marks ) {
        
        print(paste(m1, G1, m2, G2))
        
        # if ((m1 %in% c("H3K36me3", "H3K9me3", "H4K20me1") & G1 == "peak_not_TSS") | 
        #     (m2 %in% c("H3K36me3", "H3K9me3", "H4K20me1") & G2 == "peak_not_TSS"))  {
        #   
        #   pvalues <- c(pvalues, 1) } else {

        q = (nrow(all.marks.groups[all.marks.groups[, m2] == G2 & all.marks.groups[, m1] == G1, ])-1)
        m = nrow(all.marks.groups[all.marks.groups[, m2] == G2, ])
        n = nrow(all.marks.groups[all.marks.groups[, m2] != G2, ])
        k = nrow(all.marks.groups[all.marks.groups[, m1] == G1, ])
        pvalues <- c(pvalues, phyper(q, m, n, k, lower.tail = ifelse((k==0 | m==0), T, F)))
        
        # }
        
      }
      
    }
    
  }
  
  matrix.pvalues <- matrix(pvalues, ncol=45, nrow=5, byrow = T)
  matrix.fdr <- t(apply(as.data.frame(matrix.pvalues), 1, p.adjust))
  
  tmp.colnames <- c()
  for ( G3 in groups3 ) {
    
    tmp.colnames <- c(tmp.colnames, rep(G3, 9))
    
  }
  
  colnames(matrix.pvalues) <- paste(tmp.colnames, marks, sep="_")
  colnames(matrix.fdr) <- paste(tmp.colnames, marks, sep="_")
  
  rownames(matrix.pvalues) <- paste(groups1, m1, sep="_")
  rownames(matrix.fdr) <- paste(groups1, m1, sep="_")
  
  all.marks.matrix.pvalues <- rbind(all.marks.matrix.pvalues, matrix.pvalues)
  all.marks.matrix.fdr <- rbind(all.marks.matrix.fdr, matrix.fdr)
  
}


# 5. check that the matrix is triangular
all.marks.matrix.pvalues.reordered  <- all.marks.matrix.pvalues[colnames(all.marks.matrix.pvalues), ]
all.marks.matrix.fdr.reordered  <- all.marks.matrix.fdr[colnames(all.marks.matrix.fdr), ]

for ( i in 1:nrow(all.marks.matrix.pvalues.reordered)) {
  
  stopifnot(all.equal(as.numeric(all.marks.matrix.pvalues.reordered[i, ]), 
                      as.numeric(all.marks.matrix.pvalues.reordered[, i])))
  
}



# 6. plot the matrix
annotation_col_df <- data.frame(group=tmp.colnames)
rownames(annotation_col_df) <- colnames(all.marks.matrix.fdr)

annotation_row_df <- data.frame(group=rep(groups1, 9))
rownames(annotation_row_df) <- rownames(all.marks.matrix.fdr)

ann_colors = list(group = palette)

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3d.pdf", width = 14, height = 14)
pheatmap(round(-log10(all.marks.matrix.fdr.reordered+(min(all.marks.matrix.fdr.reordered[all.marks.matrix.fdr.reordered>0]))), 2),
         cluster_rows = T, cluster_cols = T,
         clustering_method = "ward.D2",
         clustering_distance_rows = "manhattan",
         clustering_distance_cols = "manhattan",
         color=c("#ddfafd", "#fde0dd", '#fcc5c0', '#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a'),
         # color=c('#fff7f3','#fde0dd','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a'),
         breaks=c(0, 1.30103, 4, 8, 10, 12, 14, 200, 500),
         legend_breaks = c(0, 1.30103, 4, 8, 10, 12, 14, 200, 500),
         legend_labels =c(0, 1.30103, 4, 8, 10, 12, 14, 200, 500),
         cellwidth = 9,
         cellheight = 9,
         border_color = NA,
         labels_row = rep(marks, 5),
         labels_col = rep(marks, 5),
         annotation_col = annotation_col_df,
         annotation_row = annotation_row_df,
         annotation_colors = ann_colors,
         annotation_names_col = F,
         annotation_names_row = F)
dev.off()
