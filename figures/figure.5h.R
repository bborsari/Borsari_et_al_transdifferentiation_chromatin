.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

#************
# LIBRARIES *
#************

library(ggplot2)
library(reshape2)
library(cowplot)


setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/ENCODE.data")

palette <- c("cluster 1" = "#F56E47", 
             "cluster 2" = "#97BF04", 
             "cluster 3"  = "#772B59")

#************
# FUNCTIONS *
#************ 

my.function1 <- function(x, y, z, class) {
  
  # x = cell line
  # y = TFs, chromatin.remodellers, RNA.pol.complex
  # z = gene.body, promoter
  
  # 1. read chromatin remodeller matrix
  m <- read.table(paste0(x, "/", y, "/", z, "/", y,
                         ".binding.plot.table.", class, ".tsv"), 
                  header = T, sep="\t")
  
  if (class == "upregulated") {
    
    # 2. retrieve fractions of genes with factor y binding in the different clusters
    m$cluster1_frac <- as.numeric(m$cluster1 / 975)
    m$cluster2_frac <- as.numeric(m$cluster2 / 1108)
    m$cluster3_frac <- as.numeric(m$cluster3 / 20)
    
    
  } else {
    
    # 2. retrieve fractions of genes with factor y binding in the different clusters
    m$cluster1_frac <- as.numeric(m$cluster1 / 2740)
    m$cluster2_frac <- as.numeric(m$cluster2 / 1232)
    m$cluster3_frac <- as.numeric(m$cluster3 / 44)
    
    
  }
  

  # 3. retrieve significant enrichments between clusters 1 and 2
  pvals <- c()
  OR <- c()
  annotation_rows <- c()
  annotation_rows2 <- c()
  
  if (class == "upregulated") {
    
    ee <- 975
    ff <- 1108
    
  } else {
    
    ee <- 2740
    ff <- 1232
    
  }
  
  for (i in 1:nrow(m)) {
    
    aa <- m[i, "cluster1"]
    bb <- ee - aa
    cc <- m[i, "cluster2"]
    dd <- ff - cc
    
    mm <- matrix(c(aa, bb, cc, dd), nrow=2, byrow = T)
    my.test <- fisher.test(mm)
    pvals <- c(pvals, my.test$p.value)
    OR <- c(OR, my.test$estimate)
    
  }
  
  m$pvals <- pvals
  m$OR <- OR
  m$fdr <- p.adjust(m$pvals, method = "fdr")
  
  for (i in 1:nrow(m)) {
    
    if (m[i, "fdr"] < 0.05 & m[i, "OR"] > 2) {
      
      annotation_rows <- c(annotation_rows, "cluster1")
      
    } else if (m[i, "fdr"] < 0.05 & m[i, "OR"] < 0.5) {
      
      annotation_rows <- c(annotation_rows, "cluster2")
      
    } else {
      
      annotation_rows <- c(annotation_rows, "ns")
      
    }
    
    
    if (m[i, "cluster3_frac"] >= 0.25) {
      
      annotation_rows2 <- c(annotation_rows2, "yes")
      
    } else {
      
      annotation_rows2 <- c(annotation_rows2, "no")
      
    }
    
  }
  
  m$annotation_rows <- as.character(annotation_rows)
  m$annotation_rows2 <- as.character(annotation_rows2)
  
  return(m)
  
}


my.function2 <- function(m, common.y) {
  
  m.common <- intersect(rownames(m), common.y)  
  m.diff <- setdiff(common.y, rownames(m))
  
  m <- m[rownames(m) %in% m.common, ]
  if (length(m.diff) > 0) {
    
    m_NA <- as.data.frame(matrix(rep(NA, ncol(m)*length(m.diff)), nrow = length(m.diff), byrow = T))
    rownames(m_NA) <- m.diff
    colnames(m_NA) <- colnames(m)
    m <- as.data.frame(rbind(as.matrix(m), as.matrix(m_NA)))
    
  }
  
  m <- m[common.y, ]
  m[, 4:6] <- apply(m[, 4:6], 2, function(x){x <- as.numeric(x)})
  return(m)
  
}


my.function3 <- function(my.list) {
  
  new.list <- list()
  n <- length(my.list)
  
  for (my.col in 1:ncol(my.list[[1]])) {
    
    tmp.df <- as.numeric(my.list[[1]][, my.col])
    
    for (my.n in 2:n) {
      
      tmp.df <- cbind(tmp.df, as.numeric(my.list[[my.n]][, my.col]))
      
    }
    
    new.list[[my.col]] <- rowMeans(tmp.df, na.rm=T)
    
  }
  
  sorted.m <- data.frame(new.list)
  return(sorted.m)
  
}


my.function4 <- function(sorted.common.y, m, cell_line){
  
  m.sub <- m[, 4:6]
  m.sub$y <- rownames(m.sub)
  m.sub.melt <- melt(m.sub)
  m.sub.melt$y <- factor(m.sub.melt$y, levels = rev(sorted.common.y))
  m.sub.melt$cell_line <- cell_line
  
  return(m.sub.melt)
  
}


my.function <- function(y, z, class, input.order = NULL) {
  
  if (y == "histone.marks") {
    
    lines <- c("GM12878", "K562", "HepG2", "MCF-7", "A549", "B.cells", "CD14.positive.monocytes")
    
  } else {
    
    lines <- c("GM12878", "K562", "HepG2", "MCF-7", "A549")
    
  }
  
  lines.list <- list()
  
  
  # 1. retrieve summary dataframes
  for (i in 1:length(lines)) {
    
    lines.list[[i]] <- my.function1(x = lines[i],
                                    y = y,
                                    z = z,
                                    class = class)
    
  }
  
  
  # 2. find factors with ChIP-seq experiments
  # in at least 2 cell lines
  all.lines <- c()
  for (i in 1:length(lines)) {
    all.lines <- c(all.lines, rownames(lines.list[[i]]))
  }
  common.y <- sort(unique(all.lines))[table(all.lines) >= 3]
  
  
  # 3. replace missing experiments within each cell line with NA
  for (i in 1:length(lines)) {
    lines.list[[i]] <- my.function2(m = lines.list[[i]],
                                    common.y = common.y)
  }

  
  # 4. check order of rownames
  for (i in 2:length(lines)) {
    stopifnot(identical(rownames(lines.list[[1]]), 
                        rownames(lines.list[[i]])))
  }
  
  
  # 5. retrieve order of rows
  lines.list2 <- list()
  for (i in 1:length(lines)) {
    
    lines.list2[[i]] <- lines.list[[i]][, 4:6]
    
  }

  
  if (z == "promoter" & is.null(input.order)) {
    
    sorted.lines <- my.function3(my.list = lines.list2)
    colnames(sorted.lines) <- colnames(lines.list[[1]][, 4:6])
    rownames(sorted.lines) <- rownames(lines.list[[1]])

    sorted.common.y <- rownames(sorted.lines[order(sorted.lines$cluster1_frac,
                                                   sorted.lines$cluster2_frac,
                                                   decreasing = T), ])
    
    
  } else {
    
    sorted.common.y <- input.order
    
  }
  
  
  
  # 6. reorder rows of dataframes according to the common y factors
  lines.sub.list <- list()
  for ( i in 1:length(lines) ) {
    lines.sub.list[[i]] <- my.function4(sorted.common.y = sorted.common.y,
                                        cell_line = lines[i],
                                        m = lines.list[[i]])
  }

  
  # 7. prepare dfs for rows annotation
  
  # 7.1. this df is to annotate significant differences in binding between cluster 1 and 2
  m_annotation_rows_all_lines <- cbind(lines.list[[1]][, "annotation_rows", drop = F],
                                       lines.list[[2]][, "annotation_rows", drop = F],
                                       lines.list[[3]][, "annotation_rows", drop = F],
                                       lines.list[[4]][, "annotation_rows", drop = F],
                                       lines.list[[5]][, "annotation_rows", drop = F])

  rownames(m_annotation_rows_all_lines) <- rownames(lines.list[[1]])
  ids <- c()
  
  for ( i in 1:nrow(m_annotation_rows_all_lines) ) {
    
    # if a gene only has an inconsistent classification across cell lines
    if (length(table(t(m_annotation_rows_all_lines[i, ]))) > 1) {
      
      ids <- c(ids, "")
      
    } else {
      
      # if a gene is significant across all cell lines
      if (names(table(t(m_annotation_rows_all_lines[i, ]))) != "ns") {
        
        ids <- c(ids, rownames(m_annotation_rows_all_lines)[i])
        
        # if a gene is not significant across all cell lines
      } else {
        
        ids <- c(ids, "")
        
      }
      
    }
    
  }
  
  m_annotation_rows_all_lines$ids <- ids
  m_annotation_rows_all_lines$alpha <- ifelse(m_annotation_rows_all_lines$ids == "", "0", "1")
  m_annotation_rows_all_lines$boldness <- ifelse(m_annotation_rows_all_lines$alpha==1, "bold", "plain")
  m_annotation_rows_all_lines$color <- ifelse(m_annotation_rows_all_lines$alpha==1, "black", "#bdbdbd")

  m_annotation_rows_all_lines <- m_annotation_rows_all_lines[sorted.common.y, ]
  
  # 7.2. this df is to annotate whether binding for cluster 3 is >= 25%
  m_annotation_rows_all_lines2 <- cbind(lines.list[[1]][, "annotation_rows2", drop = F],
                                        lines.list[[2]][, "annotation_rows2", drop = F],
                                        lines.list[[3]][, "annotation_rows2", drop = F],
                                        lines.list[[4]][, "annotation_rows2", drop = F],
                                        lines.list[[5]][, "annotation_rows2", drop = F])
  
  rownames(m_annotation_rows_all_lines2) <- rownames(lines.list[[1]])
  ids2 <- c()
  
  for ( i in 1:nrow(m_annotation_rows_all_lines2) ) {
    
    # if a gene only has an inconsistent classification across cell lines
    if (length(table(t(m_annotation_rows_all_lines2[i, ]))) > 1) {
      
      ids2 <- c(ids2, "")
      
    } else {
      
      # if a gene is significant across all cell lines
      if (names(table(t(m_annotation_rows_all_lines2[i, ]))) != "no") {
        
        ids2 <- c(ids2, "*")
        
        # if a gene is not significant across all cell lines
      } else {
        
        ids2 <- c(ids2, "")
        
      }
      
    }
    
  }
  
  
  
  
  # 8. prepare dataframe for plot
  all.sub.melt <- data.frame(stringsAsFactors = F)
  for ( i in 1:length(lines) ){
    
    tmp <- lines.sub.list[[i]]
    tmp$y <- paste0(ids2, tmp$y)
    all.sub.melt <- rbind(all.sub.melt,
                          tmp)
    
  }
  
  if ( y == "histone.marks") {
    
    all.sub.melt$cell_line <- gsub("B.cells", "B cells", all.sub.melt$cell_line)
    all.sub.melt$cell_line <- gsub("CD14.positive.monocytes", 
                                   "CD14+\nmonocytes", all.sub.melt$cell_line)
    
    all.sub.melt$cell_line <- factor(all.sub.melt$cell_line, levels =
                                       c("MCF-7", "HepG2", "A549",
                                         "GM12878", "K562", "B cells", 
                                         "CD14+\nmonocytes"))
    
  } else {
    
    
    all.sub.melt$cell_line <- factor(all.sub.melt$cell_line, levels =
                                       c("MCF-7", "HepG2", "A549",
                                         "GM12878", "K562"))
    
    
  }
  
  
  
  all.sub.melt <- merge(all.sub.melt, m_annotation_rows_all_lines[, c("ids", "alpha")], 
                        by.x = "y", by.y = "ids", all.x = T)
  all.sub.melt$alpha <- ifelse(is.na(all.sub.melt$alpha), 0, 1)
  
  old.ids <- data.frame(id1=rownames(lines.list[[1]]),
                        id2 = ids2)
  rownames(old.ids) <- old.ids$id1
  old.ids <- old.ids[sorted.common.y, ]
  sorted.common.y <- paste0(old.ids$id2, old.ids$id1)

  all.sub.melt$y <- factor(all.sub.melt$y, levels = rev(sorted.common.y))
  
  all.sub.melt$variable <- gsub("cluster1_frac", "cluster 1", all.sub.melt$variable)
  all.sub.melt$variable <- gsub("cluster2_frac", "cluster 2", all.sub.melt$variable)
  all.sub.melt$variable <- gsub("cluster3_frac", "cluster 3", all.sub.melt$variable)
  
  p <- ggplot(all.sub.melt, 
              aes(x=y, y=value*100, colour=variable)) +
    geom_point(size=2.5) +
    coord_flip() +
    facet_wrap(~cell_line, nrow = 1) +
    theme_bw() +
    theme(axis.title.x = element_text(size=15),
          axis.text.y = element_text(face = rev(m_annotation_rows_all_lines$boldness),
                                     size = 12,
                                     color = rev(m_annotation_rows_all_lines$color)),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=12, angle = 30, vjust = .5),
          strip.text.x = element_text(size=15),
          strip.background.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size=12),
          panel.border = element_rect(color="black"), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_color_manual(values = palette) +
    guides(alpha=F) +
    ylab(paste0("% of genes with peaks (", gsub("\\.", " ", z), ") - ", class, " genes")) +
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100),
                       labels = c(0, 25, 50, 75, 100)) +
    expand_limits(y=100) 

  my.legend <- get_legend(p)
  
  p <- p + guides(color = F)
  
  if (z == "promoter") {
    
    return (list(p, sorted.common.y, my.legend))
    
  } else {
    
    return(p)
    
  }
  
  
}


#********
# BEGIN *
#********


# 1. RNA.pol.complex
lop.RNA.pol.complex <- list()
RNA.pol.complex.promoter.upreg <- my.function(y="RNA.pol.complex", z="promoter", class="upregulated")
lop.RNA.pol.complex[[1]] <- RNA.pol.complex.promoter.upreg[[1]]
lop.RNA.pol.complex[[2]] <- my.function(y="RNA.pol.complex",
                            z="gene.body",
                            class="upregulated",
                            input.order = RNA.pol.complex.promoter.upreg[[2]])
lop.RNA.pol.complex[[3]] <- RNA.pol.complex.promoter.upreg[[3]]
RNA.pol.complex.promoter.downreg <- my.function(y="RNA.pol.complex", z="promoter",
                                    class="downregulated",
                                    input.order = RNA.pol.complex.promoter.upreg[[2]])
lop.RNA.pol.complex[[4]] <- RNA.pol.complex.promoter.downreg[[1]]
lop.RNA.pol.complex[[5]] <- my.function(y="RNA.pol.complex",
                            z="gene.body",
                            class="downregulated",
                            input.order = RNA.pol.complex.promoter.upreg[[2]])



pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.5h.RNA.pol.complex.promoter.pdf",
    height = 3, width = 21)
a <- plot_grid(plotlist = lop.RNA.pol.complex[c(1,4)],
               nrow = 1, ncol =2)
plot_grid(plotlist = list(a, lop.RNA.pol.complex[[3]]), nrow=1, ncol=2,
          rel_widths = c(1, 0.1))
dev.off()


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.5h.RNA.pol.complex.gb.pdf",
    height = 3, width = 21)
a <- plot_grid(plotlist = lop.RNA.pol.complex[c(2,5)],
               nrow = 1, ncol =2)
plot_grid(plotlist = list(a, lop.RNA.pol.complex[[3]]), nrow=1, ncol=2,
          rel_widths = c(1, 0.1))
dev.off()



# 2. TFs
lop.TFs <- list()
TFs.promoter.upreg <- my.function(y="TFs", z="promoter", class="upregulated")
lop.TFs[[1]] <- TFs.promoter.upreg[[1]]
lop.TFs[[2]] <- my.function(y="TFs",
                            z="gene.body",
                            class="upregulated",
                            input.order = TFs.promoter.upreg[[2]])
lop.TFs[[3]] <- TFs.promoter.upreg[[3]]
TFs.promoter.downreg <- my.function(y="TFs", z="promoter",
                                    class="downregulated",
                                    input.order = TFs.promoter.upreg[[2]])
lop.TFs[[4]] <- TFs.promoter.downreg[[1]]
lop.TFs[[5]] <- my.function(y="TFs",
                            z="gene.body",
                            class="downregulated",
                            input.order = TFs.promoter.upreg[[2]])



pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.5h.TFs.promoter.pdf",
    height = 11, width = 21)
a <- plot_grid(plotlist = lop.TFs[c(1,4)],
               nrow = 1, ncol =2)
plot_grid(plotlist = list(a, lop.TFs[[3]]), nrow=1, ncol=2,
          rel_widths = c(1, 0.1))
dev.off()


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.5h.TFs.gb.pdf",
    height = 11, width = 21)
a <- plot_grid(plotlist = lop.TFs[c(2,5)],
               nrow = 1, ncol =2)
plot_grid(plotlist = list(a, lop.TFs[[3]]), nrow=1, ncol=2,
          rel_widths = c(1, 0.1))
dev.off()




# 3. chromatin.remodellers
lop.chromatin.remodellers <- list()
chromatin.remodellers.promoter.upreg <- my.function(y="chromatin.remodellers", z="promoter", class="upregulated")
lop.chromatin.remodellers[[1]] <- chromatin.remodellers.promoter.upreg[[1]]
lop.chromatin.remodellers[[2]] <- my.function(y="chromatin.remodellers", 
                            z="gene.body", 
                            class="upregulated", 
                            input.order = chromatin.remodellers.promoter.upreg[[2]])
lop.chromatin.remodellers[[3]] <- chromatin.remodellers.promoter.upreg[[3]]
chromatin.remodellers.promoter.downreg <- my.function(y="chromatin.remodellers", z="promoter", 
                                    class="downregulated",
                                    input.order = chromatin.remodellers.promoter.upreg[[2]])
lop.chromatin.remodellers[[4]] <- chromatin.remodellers.promoter.downreg[[1]]
lop.chromatin.remodellers[[5]] <- my.function(y="chromatin.remodellers", 
                            z="gene.body", 
                            class="downregulated", 
                            input.order = chromatin.remodellers.promoter.upreg[[2]])



pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.5h.chromatin.remodellers.promoter.pdf",
    height = 5, width = 21)
a <- plot_grid(plotlist = lop.chromatin.remodellers[c(1,4)],
               nrow = 1, ncol =2)
plot_grid(plotlist = list(a, lop.chromatin.remodellers[[3]]), nrow=1, ncol=2,
          rel_widths = c(1, 0.1))
dev.off()


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.5h.chromatin.remodellers.gb.pdf",
    height = 11, width = 21)
a <- plot_grid(plotlist = lop.chromatin.remodellers[c(2,5)],
               nrow = 1, ncol =2)
plot_grid(plotlist = list(a, lop.chromatin.remodellers[[3]]), nrow=1, ncol=2,
          rel_widths = c(1, 0.1))
dev.off()





# 4. histone.marks
lop.histone.marks <- list()
histone.marks.promoter.upreg <- my.function(y="histone.marks", z="promoter", class="upregulated",
                                            input.order = c("H2AFZ", "H3K4me1", "H3K4me2",
                                                            "H3K27ac", "H3K9ac", "H3K4me3",
                                                            "H3K79me2", "H3K36me3", "H4K20me1",
                                                            "H3K9me3", "H3K27me3"))
lop.histone.marks[[1]] <- histone.marks.promoter.upreg[[1]]
lop.histone.marks[[2]] <- my.function(y="histone.marks",
                                      z="gene.body",
                                      class="upregulated",
                                      input.order = histone.marks.promoter.upreg[[2]])
lop.histone.marks[[3]] <- histone.marks.promoter.upreg[[3]]
histone.marks.promoter.downreg <- my.function(y="histone.marks", z="promoter",
                                              class="downregulated",
                                              input.order = histone.marks.promoter.upreg[[2]])
lop.histone.marks[[4]] <- histone.marks.promoter.downreg[[1]]
lop.histone.marks[[5]] <- my.function(y="histone.marks",
                                      z="gene.body",
                                      class="downregulated",
                                      input.order = histone.marks.promoter.upreg[[2]])



pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.5h.histone.marks.promoter.pdf",
    height = 8, width = 11)
a <- plot_grid(plotlist = lop.histone.marks[c(1,4)],
          nrow = 2, ncol =1)
plot_grid(plotlist = list(a, lop.histone.marks[[3]]), nrow=1, ncol=2,
          rel_widths = c(1, 0.1))
dev.off()


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.5h.histone.marks.gb.pdf",
    height = 8, width = 11)
a <- plot_grid(plotlist = lop.histone.marks[c(2,5)],
               nrow = 2, ncol =1)
plot_grid(plotlist = list(a, lop.histone.marks[[3]]), nrow=1, ncol=2,
          rel_widths = c(1, 0.1))
dev.off()

