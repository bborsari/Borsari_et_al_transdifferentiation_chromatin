.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/QC/aggregation.plots/")

#**********
# PALETTE *
#**********


palette <- c("#D3DCE0", "#9FB4C4", "#798FA6", "#4D6478", "#33475C", "black")
palette <- colorRampPalette(palette)(12)




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
  
  aggregation1 <- read.table(paste0("out/", mark, ".aggregation.plot.R1.gene.body.tsv"), h=T)
  aggregation2 <- read.table(paste0("out/", mark, ".aggregation.plot.R2.gene.body.tsv"), h=T)
  
  my.length=(nrow(aggregation1[aggregation1$sample==paste0("H000", mark, "X1.pileup_signal"), ]) - 10000)
  
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
    theme(panel.border = element_rect(color="black"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.title = element_text(size=15),
          legend.text = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.text.x = element_text(size = 15, angle = 30, vjust = .5),
          axis.title = element_text(size = 15),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          strip.text.x = element_text(size = 17),
          strip.background.x = element_blank()) +
    scale_x_continuous(breaks=c(-5000, 0, my.length, my.length+5000),
                       labels=c("-5 Kb", "0%", "100%", "+5 Kb")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = my.length, linetype = "dashed") +
    ylab("") +
    xlab("") +
    labs(color = "time-point", linetype = "replicate")
  
  if (plot.axis.title.y) {
    
    p <- p +
      ylab("pile-up signal")
    
  }
  
  if (plot.axis.title.x) {
    
    p <- p +
      xlab("distance over gene body")
    
  }
  
  
  
  return(p)
  
}


#********
# BEGIN *
#********

lop <- list()

lop[[1]] <- aggregation.plot(mark = "H3K36me3",
                             my.string = "H3K",
                             plot.axis.title.x = T,
                             plot.axis.title.y = T)

lop[[2]] <- aggregation.plot(mark = "H4K20me1",
                             my.string = "H4K",
                             plot.axis.title.x = T)



my.legend <- get_legend(lop[[1]])

lop <- lapply(lop, function(x){x <- x + guides(color=F, linetype=F)})


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.7f.pdf", 
    width = 14, height=5)
main.p <- plot_grid(plotlist = lop, nrow=1, ncol=2, align="h", scale = .9)
whole.p <- plot_grid(main.p, my.legend, nrow=1, rel_widths = c(0.9, 0.1))
whole.p
dev.off()
