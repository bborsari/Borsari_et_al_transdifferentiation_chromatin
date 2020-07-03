.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")

#************
# LIBRARIES *
#************

library(ggplot2)
library(scales)
library(ggrepel)


#************
# FUNCTIONS *
#************

my.function <- function(one.data) {
  
  one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
  one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
  one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
  one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
  one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
  one.data$frequency <- as.numeric( as.character(one.data$frequency) );
  one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
  one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
  
  return(one.data)
  
}


#********
# BEGIN *
#********


my.colnames <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability")



revigo.data <- as.data.frame(rbind(c("GO:0003002","regionalization", 0.117, 1.339, 0.362, 4.177,-4.3224,0.647,0.000),
                                   c("GO:0045893","positive regulation of transcription, DNA-templated", 0.519,-0.582, 0.549, 4.823,-4.9469,0.218,0.000),
                                   c("GO:0006357","regulation of transcription from RNA polymerase II promoter", 1.273,-0.602, 0.698, 5.213,-4.6799,0.353,0.433),
                                   c("GO:0021953","central nervous system neuron differentiation", 0.040, 1.275, 0.600, 3.706,-3.7670,0.634,0.562),
                                   c("GO:0006355","regulation of transcription, DNA-templated", 9.917,-0.595, 0.761, 6.105,-5.0521,0.226,0.622),
                                   c("GO:0009653","anatomical structure morphogenesis", 1.542, 1.423, 0.560, 5.296,-4.2757,0.640,0.656),
                                   c("GO:0007389","pattern specification process", 0.147, 1.354, 0.382, 4.274,-4.2464,0.656,0.669),
                                   c("GO:0048598","embryonic morphogenesis", 0.171, 1.364, 0.395, 4.341,-4.5834,0.644,0.679),
                                   c("GO:0060284","regulation of cell development", 0.174, 1.032, 1.023, 4.350,-3.7595,0.382,0.680)))

colnames(revigo.data) <- my.colnames
revigo.data <- my.function(one.data=revigo.data)

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.1v.pdf",
    height = 4, width = 8)
ggplot(revigo.data, aes( x=plot_X, y=plot_Y, size = -log10_p_value, label=description)) +
  geom_point(color="#737373") + 
  theme_bw() + 
  geom_text_repel(size = 4, force = 4) + 
  labs(y = "semantic space y", 
       x = "semantic space x",
       size="-log10(p-value)") + 
  coord_cartesian() +
  guides(alpha=F) +
  theme(panel.border = element_rect(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=11),
        axis.text = element_text(size=11),
        plot.title = element_blank()) 
dev.off()