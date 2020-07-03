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
revigo.data.cluster1 <- as.data.frame(rbind(c("GO:0003743","translation initiation factor activity", 0.411,-3.224,-5.375, 4.762,-4.1669,0.527,0.000),
                                            c("GO:0004842","ubiquitin-protein transferase activity", 0.352, 5.088,-1.689, 4.695,-5.3706,0.683,0.000),
                                            c("GO:0003676","nucleic acid binding",21.226,-5.098, 0.865, 6.475,-15.8697,0.534,0.149),
                                            c("GO:1901363","heterocyclic compound binding",41.115,-2.748, 4.577, 6.763,-10.7570,0.620,0.213),
                                            c("GO:0019787","ubiquitin-like protein transferase activity", 0.378, 5.107, 1.555, 4.726,-5.5817,0.683,0.219),
                                            c("GO:0003723","RNA binding", 5.283,-4.565,-2.927, 5.871,-11.2328,0.465,0.289),
                                            c("GO:0097159","organic cyclic compound binding",41.137,-4.626, 4.402, 6.763,-10.5346,0.620,0.292),
                                            c("GO:0003677","DNA binding",12.549,-4.247,-1.602, 6.247,-4.4449,0.446,0.487)))

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


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8j.cluster1.pdf", 
    width=5, height=4)
p1
dev.off()



# 2. cluster 2
revigo.data.cluster2 <- as.data.frame(rbind(c("GO:0000981","RNA polymerase II transcription factor activity, sequence-specific DNA binding", 0.640, 6.221, 0.301, 4.954,-4.0367,0.957,0.000),
                                            c("GO:0005102","receptor binding", 0.420,-6.853,-2.707, 4.772,-8.6778,0.749,0.000),
                                            c("GO:0005201","extracellular matrix structural constituent", 0.028,-0.791,-4.360, 3.601,-3.3242,0.971,0.000),
                                            c("GO:0008235","metalloexopeptidase activity", 0.229,-2.125, 3.570, 4.509,-4.2328,0.954,0.000),
                                            c("GO:0030695","GTPase regulator activity", 0.206,-3.773,-7.912, 4.463,-5.0799,0.791,0.000),
                                            c("GO:0038023","signaling receptor activity", 2.150, 1.689,-6.921, 5.481,-16.9788,0.744,0.000),
                                            c("GO:0046943","carboxylic acid transmembrane transporter activity", 0.594,-0.333, 6.163, 4.923,-3.3904,0.823,0.000),
                                            c("GO:0060089","molecular transducer activity", 2.707, 6.315,-3.568, 5.581,-14.5638,0.972,0.000),
                                            c("GO:0098772","molecular function regulator", 1.115, 1.063, 0.810, 5.196,-9.7986,0.972,0.000),
                                            c("GO:0004713","protein tyrosine kinase activity", 0.138, 4.157,-0.844, 4.288,-3.2328,0.970,0.023),
                                            c("GO:0050840","extracellular matrix binding", 0.011,-5.543, 4.191, 3.198,-3.2351,0.942,0.036),
                                            c("GO:0042277","peptide binding", 0.089, 5.304, 2.602, 4.098,-3.9469,0.920,0.042),
                                            c("GO:0005539","glycosaminoglycan binding", 0.107,-6.398, 4.589, 4.178,-3.8928,0.928,0.042),
                                            c("GO:1990837","sequence-specific double-stranded DNA binding", 0.279,-2.939, 2.126, 4.595,-3.0026,0.913,0.046),
                                            c("GO:0033218","amide binding", 0.413, 0.806,-0.895, 4.764,-4.3261,0.933,0.047),
                                            c("GO:0008289","lipid binding", 0.519, 1.631,-2.468, 4.864,-8.2013,0.932,0.048),
                                            c("GO:0005509","calcium ion binding", 0.967, 2.549, 4.008, 5.134,-7.7570,0.919,0.053),
                                            c("GO:0030246","carbohydrate binding", 0.723, 4.733,-3.159, 5.007,-5.0540,0.931,0.054),
                                            c("GO:0044877","macromolecular complex binding", 0.740, 6.159,-2.100, 5.018,-3.4330,0.931,0.054),
                                            c("GO:0005543","phospholipid binding", 0.270, 2.875, 3.041, 4.579,-6.7352,0.902,0.138),
                                            c("GO:0008201","heparin binding", 0.029,-4.694, 5.190, 3.605,-3.8210,0.919,0.194),
                                            c("GO:0043168","anion binding",20.942, 3.369, 4.235, 6.470,-3.7721,0.908,0.235),
                                            c("GO:0001848","complement binding", 0.002,-6.562,-1.065, 2.474,-3.0367,0.811,0.367),
                                            c("GO:0043236","laminin binding", 0.006,-6.227,-3.021, 2.890,-3.2381,0.802,0.391),
                                            c("GO:0043394","proteoglycan binding", 0.011,-6.721,-0.569, 3.195,-3.0888,0.788,0.409),
                                            c("GO:0019838","growth factor binding", 0.034,-6.014,-2.138, 3.677,-5.1221,0.783,0.443),
                                            c("GO:0005516","calmodulin binding", 0.072,-7.116,-1.556, 4.008,-3.0783,0.774,0.470),
                                            c("GO:0050839","cell adhesion molecule binding", 0.088,-6.622,-1.632, 4.095,-3.0835,0.771,0.477),
                                            c("GO:0038187","pattern recognition receptor activity", 0.003, 2.790,-6.529, 2.634,-3.6091,0.828,0.494),
                                            c("GO:0001228","transcriptional activator activity, RNA polymerase II transcription regulatory region sequence-specific binding", 0.079, 6.307, 0.050, 4.044,-4.3904,0.957,0.497),
                                            c("GO:0019901","protein kinase binding", 0.149,-7.114,-1.963, 4.321,-4.6676,0.748,0.498),
                                            c("GO:0030228","lipoprotein particle receptor activity", 0.004, 2.652,-6.957, 2.775,-3.4750,0.813,0.506),
                                            c("GO:0015026","coreceptor activity", 0.006, 1.565,-6.621, 2.956,-3.5560,0.809,0.521),
                                            c("GO:0003779","actin binding", 0.326,-6.880,-2.502, 4.662,-7.9747,0.746,0.533),
                                            c("GO:0048018","receptor agonist activity", 0.003,-4.413,-4.747, 2.691,-4.9393,0.612,0.541),
                                            c("GO:0042802","identical protein binding", 0.400,-6.578,-2.030, 4.751,-4.9136,0.750,0.543),
                                            c("GO:0001614","purinergic nucleotide receptor activity", 0.012, 1.049,-7.243, 3.223,-3.0367,0.783,0.546),
                                            c("GO:0016502","nucleotide receptor activity", 0.012, 1.375,-7.290, 3.223,-3.0367,0.783,0.546),
                                            c("GO:0035586","purinergic receptor activity", 0.015, 1.342,-7.048, 3.321,-3.4711,0.781,0.556),
                                            c("GO:0030545","receptor regulator activity", 0.009,-3.556,-7.834, 3.084,-5.0132,0.831,0.572),
                                            c("GO:0008092","cytoskeletal protein binding", 0.708,-6.465,-2.394, 4.999,-5.9431,0.741,0.573),
                                            c("GO:0005342","organic acid transmembrane transporter activity", 0.597, 0.119, 6.046, 4.924,-3.3904,0.844,0.576),
                                            c("GO:0015318","inorganic solute uptake transmembrane transporter activity", 0.004,-0.710, 6.168, 2.777,-3.2628,0.865,0.580),
                                            c("GO:0008238","exopeptidase activity", 0.865,-1.440, 3.125, 5.086,-3.4283,0.954,0.593),
                                            c("GO:0019899","enzyme binding", 0.618,-6.916,-2.235, 4.940,-3.7100,0.743,0.595),
                                            c("GO:0038024","cargo receptor activity", 0.046, 2.342,-6.805, 3.810,-8.2588,0.804,0.610),
                                            c("GO:0001664","G-protein coupled receptor binding", 0.066,-6.473,-3.117, 3.970,-4.5421,0.768,0.620),
                                            c("GO:0046873","metal ion transmembrane transporter activity", 0.999,-0.225, 6.030, 5.148,-3.1415,0.823,0.651),
                                            c("GO:0001540","beta-amyloid binding", 0.007, 5.174, 2.289, 3.006,-3.6162,0.925,0.660),
                                            c("GO:0015085","calcium ion transmembrane transporter activity", 0.127,-1.015, 6.161, 4.253,-3.0585,0.843,0.689),
                                            c("GO:0005088","Ras guanyl-nucleotide exchange factor activity", 0.128,-3.516,-7.997, 4.255,-3.6198,0.806,0.692)))

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

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8j.cluster2.pdf", 
    width=7, height=8)
p2
dev.off()



# 3. cluster 3
revigo.data.cluster3 <- as.data.frame(rbind(c("GO:0008009","chemokine activity", 0.016, 5.145, 2.373, 3.365,-9.0565,0.120,0.000),
                                            c("GO:0030545","receptor regulator activity", 0.009, 2.404,-5.218, 3.084,-6.8633,0.753,0.000),
                                            c("GO:0005539","glycosaminoglycan binding", 0.107,-5.337,-0.447, 4.178,-5.8697,0.746,0.033),
                                            c("GO:0042834","peptidoglycan binding", 0.065,-4.719, 3.512, 3.963,-5.5058,0.747,0.205),
                                            c("GO:0048018","receptor agonist activity", 0.003, 4.275,-0.389, 2.691,-6.0565,0.200,0.574),
                                            c("GO:0005125","cytokine activity", 0.054, 4.580, 2.983, 3.883,-7.0942,0.183,0.671)))


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

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.8j.cluster3.pdf", 
    width=5, height=4)
p3
dev.off()
