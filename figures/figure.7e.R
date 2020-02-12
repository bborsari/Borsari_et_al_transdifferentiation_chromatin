.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/QC/aggregation.plots/")

#**********
# PALETTE *
#**********


palette <- c("#662506", "#993404", "#cc4c02",
             "#ec7014", "#fe9929", "#fec44f",
             '#addd8e','#78c679','#41ab5d',
             '#238443','#006837','#004529'
)


names(palette) <- c("H000", "H003", "H006", "H009",
                    "H012", "H018", "H024", "H036",
                    "H048", "H072", "H120", "H168")

#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)


#************
# FUNCTIONS *
#************

aggregation.plot <- function(mark,
                             my.string,
                             plot.axis.title.x = F,
                             plot.axis.title.y = F) {
  
  aggregation1 <- read.table(paste0("out/", mark, ".aggregation.plot.R1.tsv"), h=T)
  aggregation2 <- read.table(paste0("out/", mark, ".aggregation.plot.R2.tsv"), h=T)
  
  aggregation1$rep <- "rep1"
  aggregation2$rep <- "rep2"
  
  aggregation.df <- rbind(aggregation1, aggregation2)
  aggregation.df$my.mark <- mark
  
  aggregation.df <- aggregation.df %>%
    separate(sample, c("tp", "other"), my.string)
  
  p <- ggplot(aggregation.df, aes(x=position, y=value, colour=tp, linetype = rep)) +
    geom_line() +
    expand_limits(y=0) +
    scale_color_manual(values = palette) +
    facet_wrap(~my.mark) +
    theme_bw() +
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          legend.position = "bottom",
          legend.title = element_text(size=15),
          legend.text = element_text(size=15),
          axis.line = element_line(colour = "black"),
          axis.text.y = element_text(size=15),
          axis.text.x = element_text(size = 15, angle = 30, vjust = .5),
          axis.title = element_text(size = 15),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          strip.text.x = element_text(size = 15)) +
    scale_x_continuous(breaks=c(-5000, -2000, 0, 2000, 5000)) +
    geom_vline(xintercept = -2000, linetype = "dashed") +
    geom_vline(xintercept = 2000, linetype = "dashed") +
    ylab("") +
    xlab("") +
    labs(color = "time-point", linetype = "replicate")
  
  if (plot.axis.title.y) {
    
    p <- p +
      ylab("pile-up signal")
    
  }
  
  if (plot.axis.title.x) {
    
    p <- p +
      xlab("distance from TSS (bp)")
    
  }
  
  
  
  return(p)
  
}


#********
# BEGIN *
#********

lop <- list()

lop[[1]] <- aggregation.plot(mark = "H3K4me3", 
                             my.string = "H3K")

lop[[2]] <- aggregation.plot(mark = "H3K9ac", 
                             my.string = "H3K")

lop[[3]] <- aggregation.plot(mark = "H3K27ac", 
                             my.string = "H3K")

lop[[4]] <- aggregation.plot(mark = "H3K4me2",
                             my.string = "H3K",
                             plot.axis.title.y = T)

lop[[5]] <- aggregation.plot(mark = "H3K9me3",
                             my.string = "H3K")

lop[[6]] <- aggregation.plot(mark = "H3K27me3",
                             my.string = "H3K")

lop[[7]] <- aggregation.plot(mark = "H3K4me1",
                             my.string = "H3K")

lop[[8]] <- aggregation.plot(mark = "H3K36me3",
                             my.string = "H3K",
                             plot.axis.title.x = T)

lop[[9]] <- aggregation.plot(mark = "H4K20me1",
                             my.string = "H4K")



my.legend <- get_legend(lop[[1]])

lop <- lapply(lop, function(x){x <- x + guides(color=F, linetype=F)})


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.7e.pdf", width = 14, height=14)
main.p <- plot_grid(plotlist = lop, nrow=3, ncol=3, align="h", scale = .9)
whole.p <- plot_grid(main.p, my.legend, nrow=2, rel_heights = c(0.9, 0.1))
whole.p
dev.off()
