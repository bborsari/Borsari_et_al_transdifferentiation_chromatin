.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks")


#**********
# PALETTE *
#**********

palette <- c("H3K9ac" = "#e7298a",
             "H3K27ac" = "#8e0152",
             "H3K4me3" = "#f1b6da",
             "H3K27me3" = "#253494",
             "H3K9me3" = "#41b6c4",
             "H3K36me3" = "#7fbc41",
             "H4K20me1" = "#276419",
             "H3K4me1" = "#ffbf00",
             "H3K4me2" = "#a67c00",
             "expression" = "#d9d9d9")


#************
# LIBRARIES *
#************

library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(cowplot)
library(reshape2)
library(dplyr)
library(plyr)
library(tidyr)
library(ggcorrplot)
library(ggpubr)


#************
# FUNCTIONS *
#************

rescale <- function(x){
  return((x -min(x))/(max(x) - min(x)))
  
}

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(my.df.melt, c("type", "variable"), .fun=summary_func,
                  "value")
  #data_sum <- rename(data_sum, c("mean" = "value"))
  return(data_sum)
}


#********
# BEGIN *
#********

# 1. read dataframes
elements <- c("H3K4me3", "H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac",
              "H3K27me3", "H3K9me3", "H3K36me3", "H4K20me1", "expression")


metadata <- read.table("H3K4me3/QN.merged/metadata.tsv", h=T, sep="\t")
upreg.genes <- rownames(metadata[metadata$final_class=="upregulation", ])


my.df <- data.frame(stringsAsFactors = F)

for (el in 1:10) {
  
  tmp <- read.table(paste0(elements[el], "/QN.merged/", elements[el], ".matrix.tsv"),
                              h=T, sep="\t")
  
  
  
  # tmp <- as.data.frame(t(apply(tmp, 1, rescale)))
  tmp$type <- elements[el]
  tmp <- tmp[rownames(tmp) %in% upreg.genes, ]
  my.df <- rbind(my.df, tmp)
  
} 

colnames(my.df) <- c("H000", "H003", "H006", "H009",
                     "H012", "H018", "H024", "H036",
                     "H048", "H072", "H120", "H168", "type")



my.df.melt <- melt(my.df)
df3 <- data_summary(my.df.melt, varname="value", 
                    groupnames=c("type", "variable"))





ggplot(my.df.melt[my.df.melt$type%in% c("H3K4me3"), ], 
       aes(x=variable, y=value, fill=type)) +
  geom_violin(fill="blue", alpha = 0.7) +
  geom_boxplot(fill="white", width = 0.25, outlier.shape = NA) +
  theme(axis.text.x = element_text(size=20, angle=45, vjust=.5),
        axis.title = element_blank(),
        axis.text.y = element_text(size=20),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_color_manual(values=palette)+
  guides(colour=F, fill=F)  +
  stat_compare_means(comparisons = list(c("H000", "H003"),
                                        c("H000", "H006"),
                                        c("H000", "H009"),
                                        c("H000", "H012"),
                                        c("H000", "H018"),
                                        c("H000", "H024"),
                                        c("H000", "H168")))
                     



# 2. rescale matrices in range 

expression.matrix <- 
H3K4me3.matrix <- read.table("H3K4me3/QN.merged/H3K4me3.matrix.tsv",
                             h=T, sep="\t")

# 2. check order of rownames
stopifnot(identical(rownames(expression.matrix),
                    rownames(H3K4me3.matrix)))

# 3. read metadata 
metadata <- read.table("H3K4me3/QN.merged/metadata.tsv", h=T, sep="\t")
metadata$gene_id <- rownames(metadata)
metadata$class <- gsub("ulation", ".", metadata$class)
metadata$class2 <- gsub("ulation", ".", metadata$class2)
metadata$class3 <- ifelse(metadata$class %in% c("downreg.", "upreg.", "peaking", "bending"),
                          as.character(metadata$class), paste0(metadata$class2, "#2"))
metadata$class4 <- ifelse(metadata$class %in% c("downreg.", "upreg.", "peaking", "bending"),
                          as.character(metadata$class), as.character(metadata$class2))
stopifnot(identical(rownames(metadata), rownames(expression.matrix)))


# 4. categories
categories.vector <- c("upreg.", "downreg.")


# 5. compute correlations

H3K4me3.sig.genes <- read.table("H3K4me3/QN.merged/H3K4me3.joint.QN.maSigPro.out.tsv",
                                h=T, sep="\t")
H3K4me3.sig.genes$gene_id <- rownames(H3K4me3.sig.genes)
H3K4me3.sig.genes <- H3K4me3.sig.genes %>% separate(gene_id, c("gene", "id"), "\\.")
rownames(H3K4me3.sig.genes) <- H3K4me3.sig.genes$gene
H3K4me3.sig.genes$gene <- NULL
H3K4me3.sig.genes$id <- NULL
H3K4me3.sig.genes <- rownames(H3K4me3.sig.genes)
H3K4me3.sig.genes <- intersect(rownames(expression.matrix), 
                               H3K4me3.sig.genes)

expression.H3K4me3.matrix <- expression.matrix[rownames(expression.matrix) %in% H3K4me3.sig.genes, ]
H3K4me3.matrix <- H3K4me3.matrix[rownames(H3K4me3.matrix) %in% H3K4me3.sig.genes, ]
stopifnot(identical(rownames(expression.H3K4me3.matrix),
                    rownames(H3K4me3.matrix)))


H3K4me3.lop <- list()

for (i in 1:length(categories.vector)) {
  
  H3K4me3.lop[[i]] <- cor.by.category(expression = expression.H3K4me3.matrix,
                                      mark = H3K4me3.matrix,
                                      metadata = metadata,
                                      category = categories.vector[i])
}



# 6. square of correlations

H3K4me3.lop2 <- list()

for (i in 1:length(H3K4me3.lop)) {
  H3K4me3.lop2[[i]] <- ggcorrplot(H3K4me3.lop[[i]],
                                  legend.title="Pearson's R",
                                  outline.col = "white",
                                  tl.col = "black",
                                  tl.cex = 0) +
    labs(title = categories.vector[i]) +
    annotate("segment", x = 1, xend = 2, y = 1, yend = 5,
             colour = "black") +
    annotate("text", x = 2, y = 5.6,
             colour = "black", label = round(diag(H3K4me3.lop[[i]])[1], 2), size = 7) +
    # annotate("segment", x = 1, xend = 2, y = 12, yend = 7.5,
    #          colour = "black") +
    # annotate("text", x = 2, y = 7,
    #          colour = "black", label = round(H3K4me3.lop[[i]][1, 12], 2), size = 7) 
    # annotate("segment", x = 12, xend = 11, y = 1, yend = 5,
    #          colour = "black") +
    # annotate("text", x = 11, y = 5.6,
    #          colour = "black", label = round(H3K4me3.lop[[i]][12, 1], 2), size = 7) 
    annotate("segment", x = 12, xend = 11, y = 12, yend = 7.5,
             colour = "black") +
    annotate("text", x = 11, y = 7,
             colour = "black", label = round(diag(H3K4me3.lop[[i]])[12], 2), size = 7)
  
  
  if (i == 1 | i==3) {
    
    H3K4me3.lop2[[i]] <- H3K4me3.lop2[[i]] + 
      scale_y_discrete(labels = c("0h", rep("", 10), "168h"))
    
  } else {
    
    H3K4me3.lop2[[i]] <- H3K4me3.lop2[[i]] + 
      scale_y_discrete(labels = rep("", 12))
    
  }
  
}

my.legend <- get_legend(H3K4me3.lop2[[1]])

H3K4me3.lop2 <- lapply(H3K4me3.lop2, function(x) {x <- x + 
  guides(fill=F) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred",
                       midpoint = 0.291065,
                       limits = c(-0.45, 0.82)) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, angle = 90, hjust = 1),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        plot.title = element_text(vjust = -0.4, size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_discrete(labels = c("0h", rep("", 10), "168h")) +
  geom_rect(size=0.75, fill=NA, colour="black",
            aes(xmin=0.5, xmax=1.5, ymin=0.5, ymax=1.5)) +
  geom_rect(size=0.75, fill=NA, colour="black",
            aes(xmin=1.5, xmax=2.5, ymin=1.5, ymax=2.5)) +
  geom_rect(size=0.75, fill=NA, colour="black",
            aes(xmin=2.5, xmax=3.5, ymin=2.5, ymax=3.5)) +
  geom_rect(size=0.75, fill=NA, colour="black",
            aes(xmin=3.5, xmax=4.5, ymin=3.5, ymax=4.5)) +
  geom_rect(size=0.75, fill=NA, colour="black",
            aes(xmin=4.5, xmax=5.5, ymin=4.5, ymax=5.5)) +
  geom_rect(size=0.75, fill=NA, colour="black",
            aes(xmin=5.5, xmax=6.5, ymin=5.5, ymax=6.5)) +
  geom_rect(size=0.75, fill=NA, colour="black",
            aes(xmin=6.5, xmax=7.5, ymin=6.5, ymax=7.5)) +
  geom_rect(size=0.75, fill=NA, colour="black",
            aes(xmin=7.5, xmax=8.5, ymin=7.5, ymax=8.5)) +
  geom_rect(size=0.75, fill=NA, colour="black",
            aes(xmin=8.5, xmax=9.5, ymin=8.5, ymax=9.5)) +
  geom_rect(size=0.75, fill=NA, colour="black",
            aes(xmin=9.5, xmax=10.5, ymin=9.5, ymax=10.5)) +
  geom_rect(size=0.75, fill=NA, colour="black",
            aes(xmin=10.5, xmax=11.5, ymin=10.5, ymax=11.5)) +
  geom_rect(size=0.75, fill=NA, colour="black",
            aes(xmin=11.5, xmax=12.5, ymin=11.5, ymax=12.5))

})


H3K4me3.plot <- plot_grid(plotlist = H3K4me3.lop2, nrow=1, ncol=2, align="h")

H3K4me3.plot <- plot_grid(H3K4me3.plot, my.legend, nrow=2, ncol=1, align = "v")

pdf("~/public_html/paper_ERC/fig.4f.pdf", height = 10)
H3K4me3.plot
dev.off()


# 7. boxplots of distributions of values

## upreg.
upreg.genes <- rownames(metadata[metadata$class4=="upreg.", ])                             
expression.H3K4me3.upreg.genes <- expression.H3K4me3.matrix[rownames(expression.H3K4me3.matrix) %in% upreg.genes, ]
H3K4me3.upreg.genes <- H3K4me3.matrix[rownames(H3K4me3.matrix) %in% upreg.genes, ]
expression.H3K4me3.upreg.genes.melt <- melt(expression.H3K4me3.upreg.genes)
expression.H3K4me3.upreg.genes.melt$type <- "expression"
H3K4me3.upreg.genes.melt <- melt(H3K4me3.upreg.genes)
H3K4me3.upreg.genes.melt$type <- "mark"
upreg.H3K4me3 <- rbind(expression.H3K4me3.upreg.genes.melt, H3K4me3.upreg.genes.melt)
upreg.H3K4me3$class <- "upreg."

my.symnum.args <- list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")) 


lop.boxplots <- list()
lop.boxplots[[1]] <- ggplot(upreg.H3K4me3[upreg.H3K4me3$type == "expression", ], 
                            aes(x=variable, y=value)) +
  geom_violin(fill="red", alpha = 0.7) +
  geom_boxplot(fill="white", width = 0.25, outlier.shape = NA) +
  labs(title="upreg. - expression") +
  ylim(0, 25) +
  ylab("log2(TPM +1)")

lop.boxplots[[2]] <- ggplot(upreg.H3K4me3[upreg.H3K4me3$type == "mark", ], 
                            aes(x=variable, y=value)) +
  geom_violin(fill="blue", alpha = 0.7) +
  geom_boxplot(fill="white", width = 0.25, outlier.shape = NA) +
  labs(title="upreg. - H3K4me3") +
  ylim(0, 17) +
  ylab("signal")

lop.boxplots <- lapply(lop.boxplots, function(x){x <- x + 
  stat_compare_means(comparisons = list(c("H000", "H003"),
                                        c("H000", "H006"),
                                        c("H000", "H009"),
                                        c("H000", "H012"),
                                        c("H000", "H018"),
                                        c("H000", "H168")),
                     method.args = list(alternative = "less"), 
                     label = "p.signif",
                     paired = T,
                     symnum.args = my.symnum.args) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size=18),
        axis.text.x = element_text(size=15, angle=30, vjust=.5),
        axis.text.y = element_text(size=15))+
  scale_x_discrete(labels = paste0(c(0,3,6,9,12,18,24,36,48,72,120,168), "h"))})

pdf("~/public_html/paper_ERC/fig.4f.bis.pdf", width=13, height=5)
plot_grid(plotlist = lop.boxplots, nrow=1, ncol=2, align="h")
dev.off()

