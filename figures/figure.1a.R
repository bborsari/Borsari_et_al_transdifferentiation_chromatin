.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


# plots of the distribution of avg and max expression for the 5 categories

#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(ggpubr)


setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/expression/QN.merged")



palette <- c("down-\nregulated" = "#2d7f89", 
             "bending" = "#7acbd5",
             "up-\nregulated" = "#89372d",
             "peaking" = "#d5847a",
             "flat" = "#8f8f8f")


#********
# BEGIN *
#********

expression.matrix <- read.table("selected.genes.rep.2.3.after.QN.merged.tsv", h=T, sep="\t")
metadata <- read.table("metadata.class2.tsv", h=T, sep="\t")
metadata$class <- gsub("regulation", "-\nregulated", metadata$class)
metadata$class2 <- gsub("regulation", "-\nregulated", metadata$class2)
metadata$class3 <- ifelse(metadata$class %in% c("down-\nregulated", "up-\nregulated", "peaking", "bending"),
                          as.character(metadata$class), paste0(metadata$class2, "#2"))
metadata$class4 <- ifelse(metadata$class %in% c("down-\nregulated", "up-\nregulated", "peaking", "bending"),
                          as.character(metadata$class), as.character(metadata$class2))
silent.genes <- read.table("../silent.genes.txt", h=F, sep="\t")
silent.genes$V1 <- gsub("\\..*", "", silent.genes$V1)
stable.genes <- setdiff(rownames(expression.matrix), rownames(metadata))
stable.genes <- setdiff(stable.genes, silent.genes$V1)
stable.genes.metadata <- data.frame(class = rep("flat", length(stable.genes)),
                                    time_point = rep(NA, length(stable.genes)),
                                    hc = rep(NA, length(stable.genes)),
                                    class2 = rep(NA, length(stable.genes)),
                                    final_class = rep(NA, length(stable.genes)),
                                    class3 = rep(NA, length(stable.genes)),
                                    class4 = rep("flat", length(stable.genes)))
rownames(stable.genes.metadata) <- stable.genes
metadata <- rbind(metadata, stable.genes.metadata)
expression.matrix <- expression.matrix[rownames(expression.matrix) %in% rownames(metadata), ]
expression.matrix <- expression.matrix[rownames(metadata), ]
stopifnot(identical(rownames(expression.matrix), rownames(metadata)))

## prepare dataframe

all.genes.exp.df <- data.frame(avg_exp = apply(expression.matrix, 1, mean),
                               max_exp = apply(expression.matrix, 1, max),
                               class = metadata$class4)

all.genes.exp.df$class <- factor(all.genes.exp.df$class, levels = c("peaking",
                                                                    "bending",
                                                                    "up-\nregulated",
                                                                    "flat",
                                                                    "down-\nregulated"))


all.genes.lop <- list()

## avg expression
all.genes.lop[[2]] <- ggplot(all.genes.exp.df, aes(x=class, y=avg_exp, fill=class)) + 
  ylab("mean log2(TPM +1)") +
  stat_compare_means(comparisons = list(c("peaking", "bending"),
                                        c("bending", "up-\nregulated"),
                                        c("up-\nregulated", "flat"),
                                        c("flat", "down-\nregulated")),
                     label = "p.format", size=5)

## max expression
all.genes.exp.df$class <- factor(all.genes.exp.df$class, levels = c("bending",
                                                                    "flat",
                                                                    "peaking",
                                                                    "down-\nregulated",
                                                                    "up-\nregulated"))

all.genes.lop[[1]] <- ggplot(all.genes.exp.df, aes(x=class, y=max_exp, fill=class)) + 
  ylab("max log2(TPM +1)") +
  stat_compare_means(comparisons = list(c("bending", "flat"),
                                        c("flat", "peaking"),
                                        c("peaking", "down-\nregulated"),
                                        c("down-\nregulated", "up-\nregulated")),
                     label = "p.format", size=5)


all.genes.lop <- lapply(all.genes.lop, function(x){x <- x + 
  geom_violin(alpha=.5, colour="black") +
  geom_boxplot(width=0.3, alpha = .7, outlier.shape = NA) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        strip.text.x = element_text(size=15),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=palette) +
  guides(fill=F)})

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1a.pdf", height = 9)
plot_grid(plotlist = all.genes.lop, nrow=2, ncol=1, align="v")
dev.off()

