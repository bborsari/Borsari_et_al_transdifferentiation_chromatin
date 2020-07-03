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

# 1. cluster 1
revigo.data.cluster1 <- as.data.frame(rbind(c("GO:0032991","macromolecular complex",14.008, 2.393,-4.237, 6.139,-11.8268,0.956,0.000),
                                            c("GO:0043226","organelle",20.788,-1.613,-4.441, 6.311,-3.7423,0.960,0.000),
                                            c("GO:0044451","nucleoplasm part", 0.895,-3.784, 5.634, 4.945,-13.8665,0.401,0.000),
                                            c("GO:1902494","catalytic complex", 3.734, 4.441, 6.098, 5.565,-16.8539,0.642,0.000),
                                            c("GO:0044424","intracellular part",35.654,-4.023, 0.155, 6.545,-13.0306,0.717,0.182),
                                            c("GO:0044441","ciliary part", 0.139,-4.291, 7.495, 4.137,-3.2790,0.564,0.379),
                                            c("GO:0005739","mitochondrion", 2.156,-4.343, 3.818, 5.327,-6.2882,0.499,0.400),
                                            c("GO:0044422","organelle part", 9.427,-5.206, 5.447, 5.967,-6.3768,0.556,0.401),
                                            c("GO:0098798","mitochondrial protein complex", 0.289, 0.162, 5.185, 4.454,-11.4921,0.329,0.416),
                                            c("GO:0098803","respiratory chain complex", 0.159, 4.904, 2.999, 4.193,-3.9393,0.681,0.418),
                                            c("GO:1990234","transferase complex", 1.223, 4.489, 5.059, 5.080,-8.8861,0.579,0.439),
                                            c("GO:0005840","ribosome", 4.198,-0.021, 4.508, 5.616,-3.6882,0.337,0.523),
                                            c("GO:0005815","microtubule organizing center", 0.350,-4.248, 6.183, 4.537,-3.1574,0.506,0.528),
                                            c("GO:0044428","nuclear part", 3.120,-3.908, 4.883, 5.487,-9.5952,0.412,0.537),
                                            c("GO:0005684","U2-type spliceosomal complex", 0.058,-0.475, 5.771, 3.759,-4.1221,0.417,0.552),
                                            c("GO:1990904","ribonucleoprotein complex", 5.291, 4.154, 6.250, 5.717,-8.7825,0.632,0.553),
                                            c("GO:0000151","ubiquitin ligase complex", 0.232, 3.250, 4.052, 4.358,-4.4622,0.506,0.604),
                                            c("GO:0043231","intracellular membrane-bounded organelle",13.760,-3.670, 4.095, 6.132,-11.5045,0.426,0.607),
                                            c("GO:0005681","spliceosomal complex", 0.250,-0.457, 4.977, 4.392,-5.3799,0.363,0.628),
                                            c("GO:0000502","proteasome complex", 0.389, 2.915, 4.409, 4.583,-3.2692,0.486,0.636),
                                            c("GO:1905368","peptidase complex", 0.452, 4.816, 4.914, 4.648,-5.2716,0.606,0.646),
                                            c("GO:0005762","mitochondrial large ribosomal subunit", 0.054,-0.128, 6.011, 3.729,-3.4168,0.340,0.654),
                                            c("GO:0043229","intracellular organelle",19.917,-3.771, 3.530, 6.292,-8.9281,0.480,0.692),
                                            c("GO:1902493","acetyltransferase complex", 0.152, 4.876, 4.317, 4.175,-3.7852,0.603,0.693),
                                            c("GO:0031461","cullin-RING ubiquitin ligase complex", 0.159, 3.381, 3.643, 4.195,-3.0888,0.517,0.695)))

colnames(revigo.data.cluster1) <- my.colnames
revigo.data.cluster1 <- my.function(one.data=revigo.data.cluster1)
ex.cluster1 <- revigo.data.cluster1[revigo.data.cluster1$dispensability < 0.3, ]

p1 <- ggplot(ex.cluster1, aes( x=plot_X, y=plot_Y, size = -log10_p_value, label=description)) +
  geom_point(color="#F56E47") + 
  theme_bw() + 
  geom_text_repel(size = 3.5) + 
  labs(y = "semantic space y", 
       x = "semantic space x",
       title = "cluster 1",
       size="-log10(p-value)") + 
  coord_cartesian() +
  guides(alpha=F) +
  theme(panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=11),
        axis.text = element_text(size=11),
        legend.position = "bottom",
        plot.title = element_text(size=15, hjust=.5)) 


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8k.cluster1.pdf", 
    width=5, height=4)
p1
dev.off()



# 2. cluster 2
revigo.data.cluster2 <- as.data.frame(rbind(c("GO:0005576","extracellular region", 2.375,-6.336,-0.147, 5.369,-14.9245,0.976,0.000),
                                            c("GO:0005623","cell",53.553, 3.249, 3.834, 6.722,-4.8239,0.989,0.000),
                                            c("GO:0005912","adherens junction", 0.146,-0.175,-5.694, 4.157,-6.2581,0.870,0.000),
                                            c("GO:0016020","membrane",61.592,-3.583, 0.276, 6.783,-17.0580,0.991,0.000),
                                            c("GO:0030054","cell junction", 0.445,-1.109, 1.177, 4.641,-7.6162,0.976,0.000),
                                            c("GO:0031982","vesicle", 1.362, 5.700,-4.666, 5.127,-15.2644,0.768,0.000),
                                            c("GO:0043235","receptor complex", 0.115, 2.215, 5.301, 4.054,-6.1169,0.963,0.000),
                                            c("GO:0044421","extracellular region part", 1.309,-0.849,-1.897, 5.110,-20.1107,0.885,0.000),
                                            c("GO:0044456","synapse part", 0.230, 0.599,-9.617, 4.355,-6.1675,0.935,0.000),
                                            c("GO:0044459","plasma membrane part", 2.405,-5.531,-4.503, 5.374,-40.3045,0.795,0.000),
                                            c("GO:0045202","synapse", 0.299,-1.216, 3.576, 4.468,-5.2472,0.976,0.000),
                                            c("GO:0045121","membrane raft", 0.076, 6.785, 0.900, 3.873,-4.8069,0.817,0.043),
                                            c("GO:0098805","whole membrane", 0.888, 1.563,-1.034, 4.941,-10.1649,0.949,0.044),
                                            c("GO:0098552","side of membrane", 0.214, 4.626, 4.471, 4.324,-7.7878,0.949,0.048),
                                            c("GO:0030496","midbody", 0.040, 5.376, 2.683, 3.590,-4.0655,0.924,0.053),
                                            c("GO:0030427","site of polarized growth", 0.091, 1.164, 2.505, 3.950,-4.9586,0.921,0.057),
                                            c("GO:0043005","neuron projection", 0.190,-0.226, 5.615, 4.271,-5.2381,0.823,0.062),
                                            c("GO:0009986","cell surface", 0.241, 3.255, 1.055, 4.375,-8.4023,0.917,0.063),
                                            c("GO:0097458","neuron part", 0.320,-1.113,-9.525, 4.498,-6.6478,0.916,0.065),
                                            c("GO:0042995","cell projection", 1.050,-2.943, 5.142, 5.014,-7.5376,0.909,0.074),
                                            c("GO:0044425","membrane part",57.394,-5.767, 2.571, 6.752,-19.7471,0.936,0.087),
                                            c("GO:0098636","protein complex involved in cell adhesion", 0.032,-3.314, 3.739, 3.501,-3.3188,0.963,0.228),
                                            c("GO:0005856","cytoskeleton", 1.542, 6.169,-4.124, 5.181,-4.5272,0.767,0.294),
                                            c("GO:0031224","intrinsic component of membrane",55.975,-6.343, 1.733, 6.741,-15.5452,0.929,0.328),
                                            c("GO:0044433","cytoplasmic vesicle part", 0.346, 4.842,-5.306, 4.532,-13.4437,0.583,0.332),
                                            c("GO:0005773","vacuole", 0.455, 4.226,-5.408, 4.651,-4.8697,0.740,0.342),
                                            c("GO:0044437","vacuolar part", 0.297, 4.197,-6.176, 4.465,-9.8268,0.649,0.383),
                                            c("GO:0005938","cell cortex", 0.252,-1.124,-6.219, 4.394,-3.6091,0.727,0.393),
                                            c("GO:0016021","integral component of membrane",55.868,-6.074, 2.055, 6.740,-13.0526,0.929,0.407),
                                            c("GO:0015629","actin cytoskeleton", 0.402, 6.441,-3.518, 4.598,-3.5800,0.788,0.481),
                                            c("GO:0005887","integral component of plasma membrane", 1.211,-4.915,-5.529, 5.076,-22.0150,0.750,0.495),
                                            c("GO:0044298","cell body membrane", 0.004,-5.508,-5.406, 2.561,-3.0367,0.824,0.498),
                                            c("GO:0001891","phagocytic cup", 0.005,-4.920,-5.005, 2.668,-3.6091,0.822,0.506),
                                            c("GO:0016324","apical plasma membrane", 0.061,-5.015,-5.399, 3.779,-4.9355,0.788,0.616),
                                            c("GO:0009897","external side of plasma membrane", 0.064,-4.597,-6.071, 3.796,-10.8508,0.795,0.618),
                                            c("GO:0001726","ruffle", 0.035, 0.029, 5.660, 3.538,-4.8356,0.843,0.634),
                                            c("GO:0005886","plasma membrane",10.510,-5.388,-4.168, 6.015,-42.7328,0.776,0.637),
                                            c("GO:0019897","extrinsic component of plasma membrane", 0.104,-4.918,-5.886, 4.009,-4.8539,0.789,0.645),
                                            c("GO:0043202","lysosomal lumen", 0.002, 3.836,-6.938, 2.326,-3.6198,0.722,0.648),
                                            c("GO:0098589","membrane region", 0.121, 6.869, 0.726, 4.077,-4.2832,0.819,0.674),
                                            c("GO:0042470","melanosome", 0.016, 4.978,-5.011, 3.188,-3.2557,0.660,0.677),
                                            c("GO:0048770","pigment granule", 0.016, 5.216,-5.541, 3.189,-3.2557,0.660,0.677),
                                            c("GO:0030667","secretory granule membrane", 0.018, 5.290,-3.276, 3.248,-7.5686,0.573,0.684),
                                            c("GO:0098590","plasma membrane region", 0.239,-4.671,-5.716, 4.372,-11.7986,0.777,0.696),
                                            c("GO:0005775","vacuolar lumen", 0.007, 3.844,-6.609, 2.846,-4.1373,0.710,0.699)))

colnames(revigo.data.cluster2) <- my.colnames
revigo.data.cluster2 <- my.function(one.data=revigo.data.cluster2)
ex.cluster2 <- revigo.data.cluster2[revigo.data.cluster2$dispensability < 0.3, ]

p2 <- ggplot(ex.cluster2, aes( x=plot_X, y=plot_Y, size = -log10_p_value, label=description)) +
  geom_point(color="#97BF04") + 
  theme_bw() + 
  geom_text_repel(size = 3) + 
  labs(y = "semantic space y", 
       x = "semantic space x",
       title = "cluster 2",
       size="-log10(p-value)") + 
  coord_cartesian() +
  guides(alpha=F) +
  theme(panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size=11),
        axis.text = element_text(size=11),
        legend.position = "bottom",
        plot.title = element_text(size=15, hjust=.5)) 

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8k.cluster2.pdf", 
    width=5, height=4)
p2
dev.off()



# 3. cluster 3
revigo.data.cluster3 <- as.data.frame(rbind(c("GO:0005576","extracellular region", 2.375, 5.013, 0.046, 5.369,-7.6925,0.455,0.000),
                                            c("GO:0005615","extracellular space", 0.538,-4.949,-0.820, 4.724,-7.5918,0.040,0.000),
                                            c("GO:0072562","blood microparticle", 0.032,-4.964, 0.730, 3.502,-5.2874,0.065,0.645)))


colnames(revigo.data.cluster3) <- my.colnames
revigo.data.cluster3 <- my.function(one.data=revigo.data.cluster3)
ex.cluster3 <- revigo.data.cluster3[revigo.data.cluster3$dispensability < 0.3, ]

p3 <- ggplot(ex.cluster3, aes( x=plot_X, y=plot_Y, size = -log10_p_value, label=description)) +
  geom_point(color="#772B59") + 
  theme_bw() + 
  geom_text_repel(size = 4) + 
  labs(y = "semantic space y", 
       x = "semantic space x",
       title = "cluster 3",
       size="-log10(p-value)") + 
  coord_cartesian() +
  guides(alpha=F) +
  theme(panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        axis.title = element_text(size=11),
        axis.text = element_text(size=11),
        plot.title = element_text(size=15, hjust=.5)) 

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8k.cluster3.pdf", 
    width=5, height=4)
p3
dev.off()
