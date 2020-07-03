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


function1 <- function(m) {
  
  titles <- c("X25_perc" = "25%",
              "X50_perc" = "50%",
              "X75_perc" = "75%",
              "X100_perc" = "100%")
  
  df <- data.frame(stringsAsFactors = F)
  
  
  for ( i in 1:8 ) {
    
    # 1. first mark
    mark1 <- lom[[i]]
    mark1$gene_id <- rownames(mark1)
    
    x <- c()
    
    for ( k in 1:8 ) {
      
      # 2. second mark
      mark2 <- lom[[k]]
      mark2$gene_id <- rownames(mark2)
      
      # 3. genes that are up-regulated and positively_correlated for the 
      # two marks
      mark12 <- merge(mark1[, c("gene_id", m)], 
                      mark2[, c("gene_id", m)], 
                      by = "gene_id")
      colnames(mark12)[2:3] <- c("mark1", "mark2")
      
      # 4. difference in time-points between 
      # up-regulation of the two marks
      mark12$diff <- mark12$mark1 - mark12$mark2
      
      type <- c()
      for ( j in 1:nrow(mark12) ) {
        
        
        if (mark12[j, "diff"] < 0) {
          
          type <- c(type, "followed by") # mark1 is followed by mark2
          
        } else if (mark12[j, "diff"] == 0) {
          
          type <- c(type, "concomitant") # mark1 is concomitant to mark2
          
        } else {
          
          type <- c(type, "anticipated by") # mark1 is anticipated by mark2
          
        }
        
      }
      
      type <- factor(type, levels = c("anticipated by", "concomitant", "followed by"))
      
      # 5. frequency of the 3 groups (anticipated by, concomitant, followed by)
      x <- c(x, round((as.vector(table(type)) / nrow(mark12))*100, 2))
      
      
      
    } 
    
    df <- rbind(df, x)
    
    
  }
  
  
  colnames(df) <- paste(rep(c("anticipated by", "concomitant", "followed by"), 8), 
                        rep(c(marks, "expression"), each=3), 
                        sep="_")
  
  # 6. reorder df to have all anticipated groups at the beginning
  rownames(df) <- c(marks, "expression")
  df <- df[marks2, c(paste(rep(c("anticipated by", "concomitant", "followed by"), each=8), 
                     rep(marks2,3), 
                     sep="_"))]
  
  
  
  # 7. prepare annotation dataframes for heatmap 
  anno_cols <- data.frame(mark1=rep(marks2,3))
  anno_rows <- data.frame(mark2=marks2)
  rownames(anno_cols) <- colnames(df)
  rownames(anno_rows) <- rownames(df)
  
  
  # 8. palette for heatmap
  ann_colors <- list(mark1 =  c("H3K4me1" = "#ffbf00",
                                "H3K4me2" = "#a67c00",
                                "H3K27ac" = "#8e0152",
                                "H3K9ac" = "#e7298a",
                                "H3K4me3" = "#f1b6da",
                                "H3K36me3" = "#7fbc41",
                                "H4K20me1" = "#276419",
                                "expression" = "gray"),
                     mark2 =  c("H3K4me1" = "#ffbf00",
                                "H3K4me2" = "#a67c00",
                                "H3K27ac" = "#8e0152",
                                "H3K9ac" = "#e7298a",
                                "H3K4me3" = "#f1b6da",
                                "H3K36me3" = "#7fbc41",
                                "H4K20me1" = "#276419",
                                "expression" = "gray"))
  
  print(df)
  
  
  pheatmap(as.matrix(df[, 1:8]),
           color=c('white', '#fef5f4', '#fef2f1', "#fde0dd", '#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a'),
           breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
           cluster_rows = F,
           cluster_cols = F,
           border_color = NA,
           cellwidth = 15,
           cellheight = 15,
           annotation_legend = F,
           annotation_names_row = F,
           annotation_names_col = F,
           show_colnames = F, 
           annotation_col = anno_cols,
           annotation_row = anno_rows,
           annotation_colors = ann_colors,
           fontsize = 15,
           main = titles[m])
  
  
  
}


#********
# BEGIN *
#********

# 1. set working directory
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


# 2. list of marks we're analyzing
marks <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3",
           "H3K36me3", "H4K20me1")

# 3. read dynamics df for each mark
lom <- list()

for ( i in 1:7 ) {
  
  # 3.1. read df
  tmp <- read.table(paste0(marks[i], "/QN.merged/ant.del.analysis/increases/", 
                           marks[i], ".25.50.75.100.tsv"),
                    h=T, sep="\t", stringsAsFactors = F)
  
  # 3.2. retrieve list of up-regulated and positively_correlated genes
  upreg.pc <- read.table(paste0(marks[i], 
                                "/QN.merged/ant.del.analysis/increases/upregulation.positively_correlated.genes.txt"),
                         h=F, sep="\t", stringsAsFactors = F)
  
  # 3.3. susbet df for the selected genes 
  lom[[i]] <- tmp[rownames(tmp) %in% upreg.pc$V1, ]
  
  
}


lom[[8]] <- read.table("H3K4me3/QN.merged/ant.del.analysis/increases/expression.25.50.75.100.tsv",
                       h=T, sep="\t", stringsAsFactors = F)



marks2 <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac",
            "expression", "H3K4me3", "H3K36me3", "H4K20me1")


# 4. make plots for 25%, 50%, 75%, 100% of upregulation
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.9r.25.pdf", width=8)
function1(m="X25_perc")
dev.off()


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.9r.50.pdf", width=8)
function1(m="X50_perc")
dev.off()


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.9r.75.pdf", width=8)
function1(m="X75_perc")
dev.off()


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.9r.100.pdf", width=8)
function1(m="X100_perc")
dev.off()