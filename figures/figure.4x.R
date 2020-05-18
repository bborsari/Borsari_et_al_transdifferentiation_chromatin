.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


#************
# LIBRARIES *
#************

library(ggplot2)
library(reshape2)
library(plyr)



#************
# FUNCTIONS *
#************

rescale <- function(x){
  return((x -min(x))/(max(x) - min(x)))
}


function1 <- function(mark) {
  
  # 1. read mark matrix
  mark.m <- read.table(paste0(mark, "/QN.merged/", mark, ".matrix.tsv"), h=T, sep="\t")
  
  # 2. rescale mark matrix in range 0-100%
  mark.m.rescaled <- t(apply(mark.m, 1, rescale))
  
  # 3. retrieve 10% least expressed genes from rescaled mark matrix
  lowest.mark.m.rescaled <- 
    mark.m.rescaled[rownames(mark.m.rescaled) %in% rownames(lowest.expression.m.rescaled), ]
  lowest.mark.m.rescaled.melt <- melt(lowest.mark.m.rescaled)
  lowest.mark.m.rescaled.melt$type <- mark
  lowest.mark.m.rescaled.melt$type2 <- "10% least expressed"

  # 4. retrieve 10% most expressed genes from rescaled mark matrix
  highest.mark.m.rescaled <- 
    mark.m.rescaled[rownames(mark.m.rescaled) %in% rownames(highest.expression.m.rescaled), ]
  highest.mark.m.rescaled.melt <- melt(highest.mark.m.rescaled)
  highest.mark.m.rescaled.melt$type <- mark
  highest.mark.m.rescaled.melt$type2 <- "10% most expressed"
  
  return(rbind(lowest.mark.m.rescaled.melt, 
               highest.mark.m.rescaled.melt))
  
}


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
  
}



#********
# BEGIN *
#********

# 1. set working directory
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


# 2. read expression matrix
expression.m <- read.table("H3K4me3/QN.merged/expression.matrix.tsv", h=T, sep="\t")


# 3. retrieve upregulated genes
metadata <- read.table("H3K4me3/QN.merged/metadata.tsv", h=T, sep="\t")
upreg.genes <- rownames(metadata[metadata$final_class == "upregulation", ])


# 4. subset expression matrix for upregulated genes
expression.m <- expression.m[rownames(expression.m) %in% upreg.genes, ]


# 5. keep only upregulated genes from cluster 2
cluster2 <- read.table("~/public_html/Borsari_et_al_transdifferentiation_chromatin/cluster.2.txt", 
                       h=F, sep="\t")
expression.m <- expression.m[rownames(expression.m) %in% cluster2$V1, ]


# 6. order expression matrix from least to most expressed genes at 0h
expression.m <- expression.m[order(expression.m$H000, decreasing = F), ]


# 7. rescale expression profiles in range 0-100%
expression.m.rescaled <- t(apply(expression.m, 1, rescale))


# 8. retrieve 
# 8.1. 10% least expressed genes
lowest.expression.m.rescaled <- expression.m.rescaled[1:111, ]
lowest.expression.m <- expression.m[1:111, ] # needed to plot distribution of expression
lowest.expression.m.rescaled.melt <- melt(lowest.expression.m.rescaled)
lowest.expression.m.rescaled.melt$type <- "expression"
lowest.expression.m.rescaled.melt$type2 <- "10% least expressed"

# 8.2. 10% most expressed genes
highest.expression.m.rescaled <- expression.m.rescaled[998:1108, ]
highest.expression.m <- expression.m[998:1108, ] # needed to plot distribution of expression
highest.expression.m.rescaled.melt <- melt(highest.expression.m.rescaled)
highest.expression.m.rescaled.melt$type <- "expression"
highest.expression.m.rescaled.melt$type2 <- "10% most expressed"

all.df <- rbind(highest.expression.m.rescaled.melt,
                lowest.expression.m.rescaled.melt)


# 8. merged df across marks
marks <- c("H3K4me1", "H3K4me2", "H3K4me3", "H3K27ac", "H3K9ac", "H3K36me3", "H4K20me1")
for ( i in marks ) {
  
  all.df <- rbind(all.df, function1(mark = i))
  
}


# 9. compute mean and SE for each mark and time-point 
all.df.summary <- summarySE(all.df, measurevar = "value", 
                            groupvars = c("Var2", "type", "type2"))


# 10. reorder marks and expression
all.df.summary$type <- factor(all.df.summary$type, levels = c("H3K4me1",
                                                              "H3K4me2",
                                                              "H3K27ac",
                                                              "expression",
                                                              "H3K9ac",
                                                              "H3K4me3",
                                                              "H3K36me3",
                                                              "H4K20me1"))

all.df.summary$type2 <- factor(all.df.summary$type2, levels = c("10% least expressed", 
                                                                "10% most expressed"))


# 10. define color palette
palette <- c("H3K9ac" = "#e7298a",
             "H3K27ac" = "#8e0152",
             "H3K4me3" = "#f1b6da",
             "H3K27me3" = "#253494",
             "H3K9me3" = "#41b6c4",
             "H3K36me3" = "#7fbc41",
             "H4K20me1" = "#276419",
             "H3K4me1" = "#ffbf00",
             "H3K4me2" = "#a67c00",
             "expression" = "black")



# 11. plot 
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.4x.pdf", 
    height = 9, width = 7.5)
ggplot(all.df.summary, aes(x=Var2, y=value*100, color=type, group=type, fill=type)) +
  geom_errorbar(aes(ymin=value*100-sd*100, ymax=value*100+sd*100), width=.1,
                position=position_dodge(0.05)) +
  geom_line() + 
  geom_point(shape = 21, color = "black", size=2.5) +  
  facet_grid(type~type2) +
  theme_bw() +
  theme(axis.title = element_text(size=11),
        axis.text.x = element_text(size=11, angle = 30, vjust = .5),
        axis.text.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15, angle = 0),
        strip.background = element_blank(),
        panel.border = element_rect(color="black"), 
        panel.grid.minor = element_blank(), 
        axis.line = element_blank()) +
  scale_x_discrete(labels = as.character(c(0,3,6,9,12,18,24,36,48,72,120,168))) +
  ylab("% of profile") +
  xlab("time (hours)") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  guides(color = F, fill=F) +
  ylim(-25, 125)
dev.off()


# 12. save list of 10% least and most expressed genes
write.table(rownames(highest.expression.m), 
            file = "~/public_html/Borsari_et_al_transdifferentiation_chromatin/genes/cluster2.highest.txt",
            row.names = F, col.names = F, quote=F, sep="\t")
write.table(rownames(lowest.expression.m), 
            file = "~/public_html/Borsari_et_al_transdifferentiation_chromatin/genes/cluster2.lowest.txt",
            row.names = F, col.names = F, quote=F, sep="\t")

