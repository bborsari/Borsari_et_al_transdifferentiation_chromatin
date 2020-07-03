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



#************************
# PALETTE & TIME-POINTS *
#************************

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


# # Pheatmap convert_annotations
# convert_annotations = function(annotation, annotation_colors){
#   new = annotation
#   for(i in 1:ncol(annotation)){
#     a = annotation[, i]
#     b = annotation_colors[[colnames(annotation)[i]]]
#     if(is.character(a) | is.factor(a)){
#       a = as.character(a)
#       
#       if(length(setdiff(setdiff(a, NA), names(b))) > 0){
#         stop(sprintf("Factor levels on variable %s do not match with annotation_colors", colnames(annotation)[i]))
#       }
#       new[, i] = b[a]
#     }
#     else{
#       a = cut(a, breaks = 100)
#       new[, i] = colorRampPalette(b)(100)[a]
#     }
#   }
#   return(as.matrix(new))
# }
# 
# # Modified pheatmap:::heatmap_motor
# heatmap_motor <- function (matrix, border_color, cellwidth, cellheight, tree_col, 
#                            tree_row, treeheight_col, treeheight_row, filename, width, 
#                            height, breaks, color, legend, annotation_row, annotation_col, 
#                            annotation_colors, annotation_legend, annotation_names_row, 
#                            annotation_names_col, main, fontsize, fontsize_row, fontsize_col, 
#                            hjust_col, vjust_col, angle_col, fmat, fontsize_number, number_color, 
#                            gaps_col, gaps_row, labels_row, labels_col, ...) 
# {
#   lo = pheatmap:::lo(coln = labels_col, rown = labels_row, nrow = nrow(matrix), 
#                      ncol = ncol(matrix), cellwidth = cellwidth, cellheight = cellheight, 
#                      treeheight_col = treeheight_col, treeheight_row = treeheight_row, 
#                      legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, 
#                      annotation_colors = annotation_colors, annotation_legend = annotation_legend, 
#                      annotation_names_row = annotation_names_row, annotation_names_col = annotation_names_col, 
#                      main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
#                      fontsize_col = fontsize_col, angle_col = angle_col, gaps_row = gaps_row, 
#                      gaps_col = gaps_col, ...)
#   res = lo$gt
#   mindim = lo$mindim
#   if (!is.na(filename)) {
#     if (is.na(height)) {
#       height = convertHeight(gtable_height(res), "inches", valueOnly = T)
#     }
#     if (is.na(width)) {
#       width = convertWidth(gtable_width(res), "inches", valueOnly = T)
#     }
#     r = regexpr("\\.[a-zA-Z]*$", filename)
#     if (r == -1) 
#       stop("Improper filename")
#     ending = substr(filename, r + 1, r + attr(r, "match.length"))
#     f = switch(ending, pdf = function(x, ...) pdf(x, ...), 
#                png = function(x, ...) png(x, units = "in", res = 300, 
#                                           ...), jpeg = function(x, ...) jpeg(x, units = "in", 
#                                                                              res = 300, ...), jpg = function(x, ...) jpeg(x, 
#                                                                                                                           units = "in", res = 300, ...), tiff = function(x, 
#                                                                                                                                                                          ...) tiff(x, units = "in", res = 300, compression = "lzw", 
#                                                                                                                                                                                    ...), bmp = function(x, ...) bmp(x, units = "in", 
#                                                                                                                                                                                                                     res = 300, ...), stop("File type should be: pdf, png, bmp, jpg, tiff"))
#     f(filename, height = height, width = width)
#     gt = heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, 
#                        border_color = border_color, tree_col = tree_col, 
#                        tree_row = tree_row, treeheight_col = treeheight_col, 
#                        treeheight_row = treeheight_row, breaks = breaks, 
#                        color = color, legend = legend, annotation_col = annotation_col, 
#                        annotation_row = annotation_row, annotation_colors = annotation_colors, 
#                        annotation_legend = annotation_legend, annotation_names_row = annotation_names_row, 
#                        annotation_names_col = annotation_names_col, filename = NA, 
#                        main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
#                        fontsize_col = fontsize_col, hjust_col = hjust_col, 
#                        vjust_col = vjust_col, angle_col = angle_col, fmat = fmat, 
#                        fontsize_number = fontsize_number, number_color = number_color, 
#                        labels_row = labels_row, labels_col = labels_col, 
#                        gaps_col = gaps_col, gaps_row = gaps_row, ...)
#     grid.draw(gt)
#     dev.off()
#     return(gt)
#   }
#   if (mindim < 3) 
#     border_color = NA
#   if (!is.na(main)) {
#     elem = pheatmap:::draw_main(main, fontsize = 1.3 * fontsize, ...)
#     res = gtable_add_grob(res, elem, t = 1, l = 3, name = "main", 
#                           clip = "off")
#   }
#   if (!pheatmap:::is.na2(tree_col) & treeheight_col != 0) {
#     elem = pheatmap:::draw_dendrogram(tree_col, gaps_col, horizontal = T)
#     res = gtable_add_grob(res, elem, t = 2, l = 3, name = "col_tree")
#   }
#   if (!pheatmap:::is.na2(tree_row) & treeheight_row != 0) {
#     elem = pheatmap:::draw_dendrogram(tree_row, gaps_row, horizontal = F)
#     res = gtable_add_grob(res, elem, t = 4, l = 1, name = "row_tree")
#   }
#   elem = pheatmap:::draw_matrix(matrix, border_color, gaps_row, gaps_col, 
#                                 fmat, fontsize_number, number_color)
#   res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", 
#                         name = "matrix")
#   if (length(labels_col) != 0) {
#     pars = list(labels_col, gaps = gaps_col, fontsize = fontsize_col, 
#                 hjust_col = hjust_col, vjust_col = vjust_col, angle_col = angle_col, 
#                 ...)
#     elem = do.call(pheatmap:::draw_colnames, pars)
#     res = gtable_add_grob(res, elem, t = 5, l = 3, clip = "off", 
#                           name = "col_names")
#   }
#   if (length(labels_row) != 0) {
#     pars = list(labels_row, gaps = gaps_row, fontsize = fontsize_row, 
#                 ...)
#     elem = do.call(pheatmap:::draw_rownames, pars)
#     res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", 
#                           name = "row_names")
#   }
#   if (!pheatmap:::is.na2(annotation_col)) {
#     converted_annotation = convert_annotations(annotation_col, 
#                                                annotation_colors)
#     elem = pheatmap:::draw_annotations(converted_annotation, border_color, 
#                                        gaps_col, fontsize, horizontal = T)
#     res = gtable_add_grob(res, elem, t = 3, l = 3, clip = "off", 
#                           name = "col_annotation")
#     if (annotation_names_col) {
#       elem = pheatmap:::draw_annotation_names(annotation_col, fontsize, 
#                                               horizontal = T)
#       res = gtable_add_grob(res, elem, t = 3, l = 4, clip = "off", 
#                             name = "col_annotation_names")
#     }
#   }
#   if (!pheatmap:::is.na2(annotation_row)) {
#     converted_annotation = convert_annotations(annotation_row, 
#                                                annotation_colors)
#     elem = pheatmap:::draw_annotations(converted_annotation, border_color, 
#                                        gaps_row, fontsize, horizontal = F)
#     res = gtable_add_grob(res, elem, t = 4, l = 2, clip = "off", 
#                           name = "row_annotation")
#     if (annotation_names_row) {
#       elem = pheatmap:::draw_annotation_names(annotation_row, fontsize, 
#                                               horizontal = F, hjust_col = hjust_col, vjust_col = vjust_col, 
#                                               angle_col = angle_col)
#       res = gtable_add_grob(res, elem, t = 5, l = 2, clip = "off", 
#                             name = "row_annotation_names")
#     }
#   }
#   annotation = c(annotation_col[length(annotation_col):1], 
#                  annotation_row[length(annotation_row):1])
#   annotation = annotation[unlist(lapply(annotation, function(x) !pheatmap:::is.na2(x)))]
#   if (length(annotation) > 0 & annotation_legend) {
#     elem = pheatmap:::draw_annotation_legend(annotation, annotation_colors, 
#                                              border_color, fontsize = fontsize, ...)
#     t = ifelse(is.null(labels_row), 4, 3)
#     res = gtable_add_grob(res, elem, t = t, l = 6, b = 5, 
#                           clip = "off", name = "annotation_legend")
#   }
#   if (!pheatmap:::is.na2(legend)) {
#     elem = pheatmap:::draw_legend(color, breaks, legend, fontsize = fontsize, 
#                                   ...)
#     t = ifelse(is.null(labels_row), 4, 3)
#     res = gtable_add_grob(res, elem, t = t, l = 5, b = 5, 
#                           clip = "off", name = "legend")
#   }
#   return(res)
# }
# 
# # Modified pheatmap:::lo    
# lo <- function (rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA, 
#                 treeheight_col, treeheight_row, legend, annotation_row, annotation_col, 
#                 annotation_colors, annotation_legend, annotation_names_row, 
#                 annotation_names_col, main, fontsize, fontsize_row, fontsize_col, 
#                 angle_col, gaps_row, gaps_col, ...) 
# {
#   if (!is.null(coln[1]) | (!pheatmap:::is.na2(annotation_row) & annotation_names_row)) {
#     if (!is.null(coln[1])) {
#       t = coln
#     }
#     else {
#       t = ""
#     }
#     tw = strwidth(t, units = "in", cex = fontsize_col/fontsize)
#     if (annotation_names_row) {
#       t = c(t, colnames(annotation_row))
#       tw = c(tw, strwidth(colnames(annotation_row), units = "in"))
#     }
#     longest_coln = which.max(tw)
#     gp = list(fontsize = ifelse(longest_coln <= length(coln), 
#                                 fontsize_col, fontsize), ...)
#     coln_height = unit(1, "grobheight", textGrob(t[longest_coln], 
#                                                  rot = angle_col, gp = do.call(gpar, gp))) + unit(10, 
#                                                                                                   "bigpts")
#   }
#   else {
#     coln_height = unit(5, "bigpts")
#   }
#   if (!is.null(rown[1])) {
#     t = rown
#     tw = strwidth(t, units = "in", cex = fontsize_row/fontsize)
#     if (annotation_names_col) {
#       t = c(t, colnames(annotation_col))
#       tw = c(tw, strwidth(colnames(annotation_col), units = "in"))
#     }
#     longest_rown = which.max(tw)
#     gp = list(fontsize = ifelse(longest_rown <= length(rown), 
#                                 fontsize_row, fontsize), ...)
#     rown_width = unit(1, "grobwidth", textGrob(t[longest_rown], 
#                                                rot = 0, gp = do.call(gpar, gp))) + unit(10, "bigpts")
#   }
#   else {
#     rown_width = unit(5, "bigpts")
#   }
#   gp = list(fontsize = fontsize, ...)
#   if (!pheatmap:::is.na2(legend)) {
#     longest_break = which.max(nchar(names(legend)))
#     longest_break = unit(1.1, "grobwidth", 
#                          textGrob(as.character(names(legend))[longest_break], 
#                                   gp = do.call(gpar, gp)))
#     title_length = unit(1.1, "grobwidth", textGrob("Scale", 
#                                                    gp = gpar(fontface = "bold", ...)))
#     legend_width = unit(12, "bigpts") + longest_break * 1.2
#     legend_width = max(title_length, legend_width)
#   }
#   else {
#     legend_width = unit(0, "bigpts")
#   }
#   if (is.na(main)) {
#     main_height = unit(0, "npc")
#   }
#   else {
#     main_height = unit(1.5, "grobheight", textGrob(main, 
#                                                    gp = gpar(fontsize = 1.3 * fontsize, ...)))
#   }
#   textheight = unit(fontsize, "bigpts")
#   if (!pheatmap:::is.na2(annotation_col)) {
#     annot_col_height = ncol(annotation_col) * (textheight + 
#                                                  unit(2, "bigpts")) + unit(2, "bigpts")
#     t = c(as.vector(as.matrix(annotation_col)), colnames(annotation_col))
#     annot_col_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], 
#                                                              gp = gpar(...))) + unit(12, "bigpts")
#     if (!annotation_legend) {
#       annot_col_legend_width = unit(0, "npc")
#     }
#   }
#   else {
#     annot_col_height = unit(0, "bigpts")
#     annot_col_legend_width = unit(0, "bigpts")
#   }
#   if (!pheatmap:::is.na2(annotation_row)) {
#     annot_row_width = ncol(annotation_row) * (textheight + 
#                                                 unit(2, "bigpts")) + unit(2, "bigpts")
#     t = c(as.vector(as.matrix(annotation_row)), colnames(annotation_row))
#     annot_row_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], 
#                                                              gp = gpar(...))) + unit(12, "bigpts")
#     if (!annotation_legend) {
#       annot_row_legend_width = unit(0, "npc")
#     }
#   }
#   else {
#     annot_row_width = unit(0, "bigpts")
#     annot_row_legend_width = unit(0, "bigpts")
#   }
#   annot_legend_width = max(annot_row_legend_width, annot_col_legend_width)
#   treeheight_col = unit(treeheight_col, "bigpts") + unit(5, 
#                                                          "bigpts")
#   treeheight_row = unit(treeheight_row, "bigpts") + unit(5, 
#                                                          "bigpts")
#   if (is.na(cellwidth)) {
#     mat_width = unit(1, "npc") - rown_width - legend_width - 
#       treeheight_row - annot_row_width - annot_legend_width
#   }
#   else {
#     mat_width = unit(cellwidth * ncol, "bigpts") + length(gaps_col) * 
#       unit(4, "bigpts")
#   }
#   if (is.na(cellheight)) {
#     mat_height = unit(1, "npc") - main_height - coln_height - 
#       treeheight_col - annot_col_height
#   }
#   else {
#     mat_height = unit(cellheight * nrow, "bigpts") + length(gaps_row) * 
#       unit(4, "bigpts")
#   }
#   gt = gtable(widths = unit.c(treeheight_row, rown_width,  
#                               mat_width, treeheight_row, legend_width, annot_legend_width), 
#               heights = unit.c(main_height, treeheight_col, annot_col_height, 
#                                mat_height, coln_height), vp = viewport(gp = do.call(gpar, 
#                                                                                     gp)))
#   cw = convertWidth(mat_width - (length(gaps_col) * unit(4, 
#                                                          "bigpts")), "bigpts", valueOnly = T)/ncol
#   ch = convertHeight(mat_height - (length(gaps_row) * unit(4, 
#                                                            "bigpts")), "bigpts", valueOnly = T)/nrow
#   mindim = min(cw, ch)
#   res = list(gt = gt, mindim = mindim)
#   return(res)
# }
# 
# # Modified pheatmap:::draw_rownames      
# draw_rownames <- function (rown, gaps, ...) 
# {
#   coord = pheatmap:::find_coordinates(length(rown), gaps)
#   y = unit(1, "npc") - (coord$coord - 0.5 * coord$size)
#   res = textGrob(rown, x = unit(-3, "bigpts"), y = y, vjust = 0.5, 
#                  hjust = 1, gp = gpar(...))
#   return(res)
# }
# 
# assignInNamespace(x="draw_rownames", value=draw_rownames, ns="pheatmap")
# assignInNamespace(x="lo", value=lo, ns="pheatmap")
# assignInNamespace(x="heatmap_motor", value=heatmap_motor, ns="pheatmap")


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
  
  f4.titles <- c("perc_25" = "25%",
                 "perc_50" = "50%",
                 "perc_75" = "75%",
                 "perc_100" = "100%")
  
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
          
          f4.aa <- nrow(f4.mark1.mark2.df[f4.mark1.mark2.df[, f4.m1] == f4.G1, ])
          f4.bb <- nrow(f4.mark1.mark2.df[f4.mark1.mark2.df[, f4.m2] == f4.G2, ])
          f4.cc <- nrow(f4.mark1.mark2.df[f4.mark1.mark2.df[, f4.m2] == f4.G2 & f4.mark1.mark2.df[, f4.m1] == f4.G1, ])
          
          f4.overlap = (f4.cc / f4.aa)
          
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
  
  
  f4.all.marks.matrix.overlap.reordered  <- f4.all.marks.matrix.overlap[colnames(f4.all.marks.matrix.overlap), ]
  
  f4.my.rows <- paste(rep(c("mark_ant", "no_diff", "exp_ant"), each=7),
                      rep(marks, 3), sep="_")
  
  f4.all.marks.matrix.overlap.reordered <- as.data.frame(f4.all.marks.matrix.overlap.reordered)[f4.my.rows, ]
  f4.all.marks.matrix.overlap.reordered <- as.data.frame(f4.all.marks.matrix.overlap.reordered)[, f4.my.rows ]
  
  f4.my.annotation.df <- data.frame(group = rep(c("mark_ant", "no_diff", "exp_ant"), each=7),
                                    mark = rep(paste(marks, "   "), 3))
  rownames(f4.my.annotation.df) <- f4.my.rows
  
  
  # write.table(f4.all.marks.matrix.overlap.reordered, 
  #             paste0("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.4v.", 
  #                    f4.degree, 
  #                    ".tsv"),
  #             row.names = T, col.names = T, sep = "\t", quote=F)
  
  f4.p <- pheatmap(as.matrix(f4.all.marks.matrix.overlap.reordered*100)[c(1:7, 15:21), ],
                   # color=c("white", '#fff7f3',"#fde0dd", '#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a'),
                   # breaks=c(0, 20, 30, 40, 50, 60, 70, 80, 90, 100),
                   color=c('white', '#fef5f4', '#fef2f1', "#fde0dd", '#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a'),
                   breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
                   cluster_rows = F,
                   cluster_cols = F,
                   border_color = NA,
                   cellwidth = 15,
                   cellheight = 15,
                   annotation_row = f4.my.annotation.df,
                   annotation_col = f4.my.annotation.df,
                   annotation_colors = list(group = palette,
                                            mark = palette2),
                   labels_row = f4.my.annotation.df[c(1:7, 15:21), "mark"],
                   # labels_col = f4.my.annotation.df[, "mark"],
                   show_colnames = F,
                   annotation_legend = F,
                   annotation_names_row = F,
                   annotation_names_col = F,
                   fontsize = 15,
                   main = f4.titles[f4.degree])
  
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
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.4v.25.pdf")
print(my.function4(f4.plot.df = my.plot.df,
                   f4.degree = "perc_25"))
dev.off()


# 3. 50% upregulation
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.4v.50.pdf")
print(my.function4(f4.plot.df = my.plot.df,
                   f4.degree = "perc_50"))
dev.off()


# 4. 75% upregulation
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.4v.75.pdf")
print(my.function4(f4.plot.df = my.plot.df,
                   f4.degree = "perc_75"))
dev.off()


# 5. 100% upregulation
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.4v.100.pdf")
print(my.function4(f4.plot.df = my.plot.df,
                   f4.degree = "perc_100"))
dev.off()
