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
revigo.data.cluster1 <- as.data.frame(rbind(c("GO:0006396","RNA processing", 3.210,-5.638, 0.964, 5.615,-11.2396,0.595,0.000),
                                            c("GO:0008152","metabolic process",75.387, 4.840,-0.287, 6.986,-10.0315,0.993,0.000),
                                            c("GO:0006325","chromatin organization", 0.668, 2.824, 5.612, 4.933,-6.4750,0.797,0.039),
                                            c("GO:0006457","protein folding", 0.903, 4.210,-1.974, 5.064,-3.9586,0.912,0.040),
                                            c("GO:0022900","electron transport chain", 0.564, 4.791, 1.703, 4.860,-4.0570,0.849,0.069),
                                            c("GO:0044237","cellular metabolic process",53.061,-1.796,-7.175, 6.833,-12.6737,0.847,0.080),
                                            c("GO:0006807","nitrogen compound metabolic process",38.744, 5.378,-1.050, 6.696,-6.9626,0.925,0.088),
                                            c("GO:1901360","organic cyclic compound metabolic process",30.324, 0.662,-5.697, 6.590,-10.7595,0.844,0.097),
                                            c("GO:0071704","organic substance metabolic process",58.357, 3.763,-4.612, 6.874,-7.8827,0.921,0.119),
                                            c("GO:0044238","primary metabolic process",53.743, 3.570,-3.951, 6.839,-7.6234,0.922,0.120),
                                            c("GO:0070647","protein modification by small protein conjugation or removal", 0.821,-0.953, 2.505, 5.023,-8.9431,0.682,0.171),
                                            c("GO:0006725","cellular aromatic compound metabolic process",29.628,-3.969,-6.461, 6.580,-12.2857,0.792,0.176),
                                            c("GO:0044265","cellular macromolecule catabolic process", 1.268,-1.009, 0.208, 5.211,-4.5229,0.652,0.180),
                                            c("GO:0043412","macromolecule modification", 9.785,-1.870,-2.500, 6.099,-3.9788,0.748,0.195),
                                            c("GO:0043170","macromolecule metabolic process",39.491,-0.050,-5.636, 6.705,-8.4711,0.837,0.211),
                                            c("GO:0043487","regulation of RNA stability", 0.021,-4.415, 6.994, 3.421,-3.6716,0.774,0.219),
                                            c("GO:0046483","heterocycle metabolic process",29.664,-4.375,-6.528, 6.580,-12.2557,0.791,0.245),
                                            c("GO:0006368","transcription elongation from RNA polymerase II promoter", 0.082,-6.359, 2.959, 4.021,-3.7258,0.683,0.293),
                                            c("GO:0010608","posttranscriptional regulation of gene expression", 0.719,-4.971, 5.199, 4.965,-3.8894,0.707,0.312),
                                            c("GO:0006354","DNA-templated transcription, elongation", 0.202,-5.984, 1.906, 4.413,-3.8729,0.667,0.320),
                                            c("GO:0006281","DNA repair", 2.234,-7.280,-0.084, 5.457,-6.0565,0.601,0.334),
                                            c("GO:0006414","translational elongation", 0.777,-3.395, 2.415, 4.999,-4.5560,0.623,0.345),
                                            c("GO:0006397","mRNA processing", 0.561,-6.149, 2.494, 4.857,-11.4437,0.585,0.358),
                                            c("GO:0006352","DNA-templated transcription, initiation", 0.766,-5.741, 1.400, 4.992,-3.6383,0.628,0.371),
                                            c("GO:0016071","mRNA metabolic process", 0.798,-7.674, 0.556, 5.010,-10.3675,0.662,0.373),
                                            c("GO:0070125","mitochondrial translational elongation", 0.009,-0.830, 4.990, 3.073,-4.3161,0.662,0.382),
                                            c("GO:0090304","nucleic acid metabolic process",21.449,-5.332,-1.490, 6.440,-15.3125,0.592,0.391),
                                            c("GO:0034645","cellular macromolecule biosynthetic process",19.291,-4.039,-0.924, 6.394,-4.4056,0.596,0.404),
                                            c("GO:0044267","cellular protein metabolic process",14.293,-2.984,-0.883, 6.263,-3.8239,0.609,0.408),
                                            c("GO:0044260","cellular macromolecule metabolic process",34.276,-3.950,-2.914, 6.643,-7.0092,0.659,0.431),
                                            c("GO:0006399","tRNA metabolic process", 2.495,-6.676,-0.100, 5.505,-3.6478,0.622,0.433),
                                            c("GO:0034641","cellular nitrogen compound metabolic process",34.137,-6.472,-3.622, 6.641,-11.4559,0.709,0.444),
                                            c("GO:0034660","ncRNA metabolic process", 3.407,-6.275,-0.063, 5.641,-4.8633,0.609,0.453),
                                            c("GO:0016567","protein ubiquitination", 0.523,-0.465, 2.383, 4.827,-7.9245,0.655,0.470),
                                            c("GO:0009059","macromolecule biosynthetic process",19.548,-3.223,-0.991, 6.399,-3.5952,0.686,0.506),
                                            c("GO:0019219","regulation of nucleobase-containing compound metabolic process",10.258,-7.490,-0.550, 6.119,-4.9747,0.617,0.519),
                                            c("GO:0044249","cellular biosynthetic process",30.048,-4.695,-4.534, 6.586,-3.5670,0.730,0.556),
                                            c("GO:0043933","macromolecular complex subunit organization", 2.371, 2.563, 6.020, 5.483,-5.7212,0.811,0.557),
                                            c("GO:0043618","regulation of transcription from RNA polymerase II promoter in response to stress", 0.026,-7.311, 3.174, 3.526,-3.6038,0.651,0.567),
                                            c("GO:0031123","RNA 3'-end processing", 0.145,-6.172, 3.463, 4.271,-3.6308,0.654,0.568),
                                            c("GO:0043620","regulation of DNA-templated transcription in response to stress", 0.027,-7.640, 2.833, 3.543,-4.1221,0.657,0.568),
                                            c("GO:0016070","RNA metabolic process",15.951,-5.394,-0.861, 6.311,-12.6990,0.562,0.577),
                                            c("GO:0006139","nucleobase-containing compound metabolic process",26.547,-6.359,-2.376, 6.532,-13.4425,0.649,0.597),
                                            c("GO:0031146","SCF-dependent proteasomal ubiquitin-dependent protein catabolic process", 0.022, 0.856, 3.196, 3.447,-4.0168,0.680,0.603),
                                            c("GO:0008380","RNA splicing", 0.413,-6.651, 2.187, 4.725,-8.4868,0.626,0.624),
                                            c("GO:0006283","transcription-coupled nucleotide-excision repair", 0.038,-9.431, 0.518, 3.690,-4.0168,0.712,0.627),
                                            c("GO:0006338","chromatin remodeling", 0.137, 2.998, 5.443, 4.244,-3.6778,0.817,0.630),
                                            c("GO:0006415","translational termination", 0.201,-1.682, 3.619, 4.411,-3.9547,0.572,0.651),
                                            c("GO:0033108","mitochondrial respiratory chain complex assembly", 0.061, 2.621, 5.650, 3.896,-4.1979,0.808,0.654),
                                            c("GO:0061013","regulation of mRNA catabolic process", 0.010,-6.399, 4.854, 3.096,-4.0899,0.620,0.665),
                                            c("GO:0032774","RNA biosynthetic process",10.925,-5.243,-0.162, 6.147,-4.0205,0.517,0.680)))

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


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3x.cluster1.pdf", 
    width=7, height=5)
p1
dev.off()



# 2. cluster 2
revigo.data.cluster2 <- as.data.frame(rbind(c("GO:0001775","cell activation", 0.171, 0.965, 2.806, 4.341,-16.3458,0.884,0.000),
                                            c("GO:0002376","immune system process", 0.600, 3.211,-0.029, 4.886,-19.1068,0.994,0.000),
                                            c("GO:0006023","aminoglycan biosynthetic process", 0.576, 0.544,-0.673, 4.869,-4.3468,0.973,0.000),
                                            c("GO:0007155","cell adhesion", 0.544, 2.155,-0.841, 4.844,-17.7545,0.927,0.000),
                                            c("GO:0022610","biological adhesion", 0.550, 3.607,-3.418, 4.849,-18.4776,0.994,0.000),
                                            c("GO:0023052","signaling", 6.765, 4.464,-0.989, 5.939,-5.8447,0.995,0.000),
                                            c("GO:0032501","multicellular organismal process", 2.373, 2.789,-2.525, 5.483,-13.5452,0.994,0.000),
                                            c("GO:0032502","developmental process", 2.812, 2.526,-1.485, 5.557,-13.5186,0.994,0.000),
                                            c("GO:0032879","regulation of localization", 0.726,-4.356, 3.421, 4.969,-19.1979,0.682,0.000),
                                            c("GO:0040011","locomotion", 0.997, 1.420,-0.172, 5.107,-16.2240,0.994,0.000),
                                            c("GO:0050896","response to stimulus",12.210, 4.042,-0.428, 6.195,-7.6216,0.995,0.000),
                                            c("GO:0051179","localization",18.495, 1.027,-0.351, 6.375,-5.8761,0.995,0.000),
                                            c("GO:0065007","biological regulation",20.498, 2.913,-0.303, 6.420,-10.7773,0.995,0.000),
                                            c("GO:0007154","cell communication", 7.219, 1.753,-2.000, 5.967,-6.9626,0.959,0.036),
                                            c("GO:0007163","establishment or maintenance of cell polarity", 0.099, 1.589, 2.850, 4.102,-4.1475,0.888,0.141),
                                            c("GO:0045321","leukocyte activation", 0.145, 3.097, 2.066, 4.269,-14.7190,0.701,0.145),
                                            c("GO:0006928","movement of cell or subcellular component", 0.973, 0.893, 3.290, 5.097,-12.6073,0.868,0.170),
                                            c("GO:0000910","cytokinesis", 0.315, 1.978, 3.290, 4.606,-4.3270,0.873,0.180),
                                            c("GO:0030029","actin filament-based process", 0.398, 0.889, 3.062, 4.708,-4.6946,0.877,0.184),
                                            c("GO:0098657","import into cell", 0.015, 3.885, 3.646, 3.288,-9.2291,0.940,0.200),
                                            c("GO:1902903","regulation of supramolecular fiber organization", 0.166,-4.839, 2.350, 4.329,-4.2652,0.697,0.234),
                                            c("GO:0032970","regulation of actin filament-based process", 0.179,-5.311, 2.432, 4.362,-5.4559,0.701,0.235),
                                            c("GO:0006826","iron ion transport", 0.133, 4.424, 3.930, 4.233,-5.4449,0.916,0.238),
                                            c("GO:0042127","regulation of cell proliferation", 0.313,-7.203, 0.328, 4.603,-11.4737,0.675,0.248),
                                            c("GO:0048878","chemical homeostasis", 0.543,-7.099, 1.570, 4.843,-6.5544,0.680,0.248),
                                            c("GO:0010941","regulation of cell death", 0.344,-5.277, 2.196, 4.645,-7.3233,0.667,0.250),
                                            c("GO:0006897","endocytosis", 0.235, 4.317, 4.218, 4.480,-7.9914,0.902,0.251),
                                            c("GO:0001932","regulation of protein phosphorylation", 0.430,-7.720,-0.727, 4.742,-7.8962,0.630,0.255),
                                            c("GO:0051239","regulation of multicellular organismal process", 0.628,-7.922, 0.922, 4.906,-15.8210,0.637,0.265),
                                            c("GO:0051049","regulation of transport", 0.529,-3.983, 3.299, 4.832,-12.7352,0.616,0.271),
                                            c("GO:0048870","cell motility", 0.633, 2.037, 4.982, 4.910,-14.8996,0.754,0.275),
                                            c("GO:2000377","regulation of reactive oxygen species metabolic process", 0.031,-6.222, 0.455, 3.595,-3.7447,0.745,0.280),
                                            c("GO:0048583","regulation of response to stimulus", 1.120,-5.596,-3.516, 5.158,-18.4260,0.617,0.281),
                                            c("GO:0072376","protein activation cascade", 0.018,-1.702,-6.739, 3.367,-3.4750,0.866,0.282),
                                            c("GO:0032940","secretion by cell", 0.763, 2.317, 4.926, 4.991,-11.8928,0.791,0.283),
                                            c("GO:0010646","regulation of cell \ncommunication", 0.929,-7.296,-0.829, 5.076,-16.3605,0.667,0.289),
                                            c("GO:0023051","regulation of signaling", 0.934,-7.694,-0.944, 5.079,-16.6556,0.673,0.289),
                                            c("GO:0050793","regulation of developmental process", 1.205,-7.325, 1.051, 5.189,-14.0655,0.623,0.297),
                                            c("GO:0016192","vesicle-mediated transport", 1.085, 3.856, 4.458, 5.144,-7.9318,0.919,0.304),
                                            c("GO:0099024","plasma membrane invagination", 0.015,-1.025, 2.257, 3.293,-3.6819,0.956,0.321),
                                            c("GO:0065008","regulation of biological quality", 3.395,-7.420, 0.403, 5.639,-11.3575,0.683,0.323),
                                            c("GO:0050790","regulation of catalytic activity", 1.575,-7.701, 0.579, 5.306,-9.7825,0.667,0.334),
                                            c("GO:0065009","regulation of molecular function", 1.726,-7.398, 0.305, 5.345,-9.5969,0.702,0.339),
                                            c("GO:0048518","positive regulation of biological process", 1.744,-8.228, 0.521, 5.350,-9.7773,0.690,0.339),
                                            c("GO:0048519","negative regulation of biological process", 1.984,-8.260, 0.469, 5.406,-7.3002,0.687,0.345),
                                            c("GO:0051174","regulation of phosphorus metabolic process", 0.580,-6.844, 0.298, 4.872,-6.8570,0.693,0.358),
                                            c("GO:0043207","response to external biotic stimulus", 0.300,-0.888,-7.354, 4.585,-8.1337,0.817,0.358),
                                            c("GO:0009607","response to biotic stimulus", 0.342,-0.799,-6.971, 4.643,-7.7670,0.859,0.363),
                                            c("GO:0048584","positive regulation of response to stimulus", 0.461,-5.414,-3.878, 4.772,-13.4711,0.525,0.374),
                                            c("GO:0009719","response to endogenous stimulus", 0.526,-1.017,-6.807, 4.829,-4.6055,0.854,0.379),
                                            c("GO:0006952","defense response", 0.568,-0.898,-6.731, 4.863,-9.2725,0.846,0.382),
                                            c("GO:0007259","JAK-STAT cascade", 0.033,-4.746,-4.368, 3.630,-4.7471,0.665,0.396),
                                            c("GO:0097696","STAT cascade", 0.033,-4.587,-3.873, 3.631,-4.7471,0.665,0.396),
                                            c("GO:0051246","regulation of protein metabolic process", 1.551,-7.742,-0.462, 5.299,-4.9508,0.674,0.402),
                                            c("GO:0070887","cellular response to chemical stimulus", 1.007,-1.237,-7.112, 5.111,-8.2381,0.792,0.406),
                                            c("GO:0007187","G-protein coupled receptor signaling pathway, coupled to cyclic nucleotide second messenger", 0.046,-4.852,-4.008, 3.771,-6.2007,0.660,0.406),
                                            c("GO:0050804","modulation of synaptic transmission", 0.057,-5.103, 1.085, 3.866,-7.9245,0.677,0.414),
                                            c("GO:1904645","response to beta-amyloid", 0.001, 0.222,-6.536, 2.057,-4.6757,0.875,0.416),
                                            c("GO:0009605","response to external stimulus", 1.370,-1.024,-7.098, 5.245,-9.7144,0.842,0.420),
                                            c("GO:0019932","second-messenger-mediated signaling", 0.079,-5.064,-4.062, 4.005,-5.6345,0.646,0.425),
                                            c("GO:0072512","trivalent inorganic cation transport", 0.014, 4.135, 3.622, 3.258,-3.9957,0.930,0.429),
                                            c("GO:0051234","establishment of localization",17.756, 4.119, 4.581, 6.358,-7.3904,0.897,0.440),
                                            c("GO:0051128","regulation of cellular component organization", 1.586,-6.709, 0.789, 5.308,-3.6440,0.667,0.462),
                                            c("GO:0042221","response to chemical", 3.071,-1.129,-6.962, 5.595,-7.2950,0.830,0.475),
                                            c("GO:0007167","enzyme linked receptor protein signaling pathway", 0.279,-5.075,-3.554, 4.554,-7.0227,0.618,0.476),
                                            c("GO:0001508","action potential", 0.024,-6.756, 2.327, 3.485,-3.4711,0.750,0.482),
                                            c("GO:0006954","inflammatory response", 0.110,-0.613,-6.744, 4.151,-8.9066,0.854,0.491),
                                            c("GO:0007267","cell-cell signaling", 0.407,-2.047,-1.659, 4.718,-4.7799,0.814,0.493),
                                            c("GO:0009611","response to wounding", 0.127,-0.474,-6.984, 4.212,-4.3307,0.862,0.497),
                                            c("GO:0050789","regulation of biological process",19.373,-7.478, 0.320, 6.395,-9.7932,0.623,0.502),
                                            c("GO:0007264","small GTPase mediated signal transduction", 0.485,-5.206,-3.567, 4.794,-4.0438,0.599,0.502),
                                            c("GO:0001505","regulation of neurotransmitter levels", 0.055,-7.047, 2.176, 3.848,-5.8827,0.737,0.514),
                                            c("GO:0050878","regulation of body fluid levels", 0.074,-7.510, 2.264, 3.976,-6.4828,0.732,0.526),
                                            c("GO:2001023","regulation of response to drug", 0.001,-3.588,-5.187, 2.155,-5.3270,0.704,0.529),
                                            c("GO:0007166","cell surface receptor signaling pathway", 0.920,-5.217,-3.221, 5.072,-15.6478,0.585,0.536),
                                            c("GO:0007186","G-protein coupled receptor signaling pathway", 0.882,-5.109,-3.203, 5.054,-12.8210,0.586,0.538),
                                            c("GO:0031099","regeneration", 0.031,-0.339, 0.728, 3.600,-3.7496,0.863,0.555),
                                            c("GO:0046903","secretion", 0.810, 3.823, 4.559, 5.017,-12.6517,0.852,0.557),
                                            c("GO:0032535","regulation of cellular component size", 0.179,-5.277, 2.758, 4.362,-4.0820,0.675,0.566),
                                            c("GO:0090066","regulation of anatomical structure size", 0.216,-7.696, 2.122, 4.444,-6.0846,0.711,0.576),
                                            c("GO:1904646","cellular response to beta-amyloid", 0.000, 0.229,-6.400, 1.301,-4.1180,0.868,0.586),
                                            c("GO:0006811","ion transport", 5.344, 4.119, 4.638, 5.836,-4.3862,0.907,0.591),
                                            c("GO:0007188","adenylate cyclase-modulating G-protein coupled receptor signaling pathway", 0.043,-4.779,-4.014, 3.743,-4.0737,0.660,0.608),
                                            c("GO:0071495","cellular response to endogenous stimulus", 0.402,-0.550,-6.540, 4.712,-3.9872,0.850,0.609),
                                            c("GO:0030001","metal ion transport", 1.677, 3.694, 4.360, 5.333,-4.5229,0.903,0.620),
                                            c("GO:0035556","intracellular signal transduction", 4.000,-5.296,-3.025, 5.710,-6.9393,0.531,0.641),
                                            c("GO:0042493","response to drug", 0.266,-0.677,-7.295, 4.534,-3.5817,0.837,0.641),
                                            c("GO:0051241","negative regulation of multicellular organismal process", 0.211,-7.588, 1.055, 4.432,-14.5086,0.579,0.650),
                                            c("GO:0090183","regulation of kidney development", 0.010,-5.770, 1.226, 3.090,-4.3401,0.660,0.654),
                                            c("GO:1901700","response to oxygen-containing compound", 0.503,-0.927,-7.378, 4.810,-7.4271,0.829,0.683),
                                            c("GO:0050818","regulation of coagulation", 0.020,-6.622, 1.083, 3.402,-3.9355,0.680,0.686),
                                            c("GO:0051952","regulation of amine transport", 0.012,-3.554, 4.281, 3.199,-4.1500,0.695,0.686),
                                            c("GO:0002819","regulation of adaptive immune response", 0.025,-3.969,-4.143, 3.513,-3.7190,0.553,0.686),
                                            c("GO:0006898","receptor-mediated endocytosis", 0.095, 4.364, 4.136, 4.086,-5.6126,0.905,0.688),
                                            c("GO:1904705","regulation of vascular smooth muscle cell proliferation", 0.004,-4.997, 0.721, 2.665,-3.6556,0.750,0.689),
                                            c("GO:0043299","leukocyte degranulation", 0.014, 2.889, 4.105, 3.252,-8.0915,0.697,0.695),
                                            c("GO:2000379","positive regulation of reactive oxygen species metabolic process", 0.016,-6.103, 0.133, 3.301,-3.8210,0.679,0.695),
                                            c("GO:0061640","cytoskeleton-dependent cytokinesis", 0.080, 1.700, 2.590, 4.014,-4.3251,0.883,0.699)))

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

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3x.cluster2.pdf", 
    width=7, height=8)
p2
dev.off()

                                      

# 3. cluster 3
revigo.data.cluster3 <- as.data.frame(rbind(c("GO:0070098","chemokine-mediated signaling pathway", 0.018, 1.118, 2.856, 3.375,-8.4522,0.446,0.000),
                                            c("GO:0061037","negative regulation of cartilage development", 0.005,-3.774,-3.561, 2.830,-4.0458,0.872,0.143),
                                            c("GO:0019730","antimicrobial humoral response", 0.010, 5.445, 3.396, 3.128,-7.0809,0.203,0.208),
                                            c("GO:0006952","defense response", 0.568, 4.951, 0.265, 4.863,-5.9031,0.524,0.268),
                                            c("GO:0070374","positive regulation of ERK1 and ERK2 cascade", 0.034, 0.706, 1.307, 3.636,-4.0264,0.523,0.306),
                                            c("GO:0009607","response to biotic stimulus", 0.342, 3.652,-0.627, 4.643,-3.9208,0.573,0.340),
                                            c("GO:0007267","cell-cell signaling", 0.407,-3.195, 5.061, 4.718,-3.9706,0.786,0.372),
                                            c("GO:0006954","inflammatory response", 0.110, 5.390,-0.716, 4.151,-6.3420,0.516,0.491),
                                            c("GO:0030593","neutrophil chemotaxis", 0.015, 3.693, 3.345, 3.272,-6.7620,0.134,0.565),
                                            c("GO:0006959","humoral immune response", 0.035, 5.323, 2.661, 3.650,-6.2396,0.312,0.663)))


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

pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/single_figures/fig.3x.cluster3.pdf", 
    width=5, height=4)
p3
dev.off()
