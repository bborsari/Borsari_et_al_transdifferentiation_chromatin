.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/QC/correlations/out")

#************
# LIBRARIES *
#************

library(ggplot2)
library(ggrepel)


#**********
# PALETTE *
#**********

palette <- c("H3K9ac" = "#a50f15",
             "H3K27ac" = "#ff7b7b",
             "H3K4me3" = "#f2629c",
             "H3K27me3" = "#253494",
             "H3K9me3" = "#6baed6",
             "H3K36me3" = "#0f9200",
             "H4K20me1" = "#4ae54a",
             "H3K4me1" = "#ffbf00",
             "H3K4me2" = "#a67c00")


#********
# BEGIN *
#********

# 1. import dataframes of cc
df <- data.frame()

for (i in 1:9) {
  
  input.df <- read.table(paste0(names(palette)[i], ".cc.tsv"), h=T) 
  df <- rbind(df, input.df)
  
}

# 2. keep Pearson's cc for 500 bp window only
df <- df[df$window==500 & df$cc == "pearson", ]

# 3. substiute label for time-point
df$time <- gsub("H000", "0", df$time)
df$time <- gsub("H003", "3", df$time)
df$time <- gsub("H006", "6", df$time)
df$time <- gsub("H009", "9", df$time)
df$time <- gsub("H012", "12", df$time)
df$time <- gsub("H018", "18", df$time)
df$time <- gsub("H024", "24", df$time)
df$time <- gsub("H036", "36", df$time)
df$time <- gsub("H048", "48", df$time)
df$time <- gsub("H072", "72", df$time)
df$time <- gsub("H120", "120", df$time)
df$time <- gsub("H168", "168", df$time)

# 4. plot
df$pos <- 1
pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.7c.pdf", width = 11, height = 3)
set.seed(123)
ggplot(df,
       aes(x=pos, y=value,
           color=mark))  +
  geom_boxplot(position = position_dodge(width=0.9), alpha=.5, width = 0.3, color = "gray")+
  ylab("Pearson's R") +
  theme_bw() +
  facet_grid(~mark) +
  expand_limits(y=1) +
  scale_color_manual(values = palette) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.text.y = element_text(size=15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size=15)) +
  guides(colour=F) +
  geom_text_repel(aes(label=time),size=4, 
                  position = position_jitterdodge(dodge.width = 1),
                  segment.alpha = 0.5,
                  force = 1.5)
dev.off()
