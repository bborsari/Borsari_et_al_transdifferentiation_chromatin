#!/usr/bin/env Rscript

.libPaths("/nfs/users2/rg/bborsari/software/R-3.4.1/library")



#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option( c ("-u", "--universe"), 
               help = "A list of gene identifiers to define the universe (ensEMBL ids), w/o header. Can be 'stdin'." ),
  
  make_option( c ("-G", "--genes"), default = 'stdin',
               help = "A list of gene identifiers for the foreground (ensEMBL ids), w/o header [default = %default]." ),
  
  make_option( c ("-c", "--category"), default = "BP",
               help = "Choose the GO category < BP | MF | CC > [default = %default]." ),
  
  make_option( c ("-s", "--species"), default = "human", 
               help = "Choose the species < human | mouse | dmel > [default = %default]." ),
  
  make_option( c ("-o", "--output"), default = "out",
               help = "Additional tags for output [default = %default]." ),
  
  make_option( c ("-d", "--directory"), default="./", 
               help = "Directory for the output [default = %default]."),
  
  make_option( c ("-p", "--pseudocount"), type = "numeric", default = 1e-300,
               help = "The pseudocount to use for the barplot [default = %default]." ),
  
  make_option( c ("-f", "--fdr_cutoff"), type = "numeric", default = 0.05,
               help = "The fdr cutoff to apply [default = %default]." ),
  
  make_option( c ("-t", "--title"), default = NULL,
               help = "The title for the barplot [default = %default]."),
               
  make_option( c("--textsize"), default = 15, type = 'numeric',
               help = "Size of plot text [default = %default]." )

)


parser <- OptionParser(
  usage = "%prog [options] file", 
  option_list=option_list,
  description = "\nProvided a universe and a list of gene ids, it returns a GO enrichment analysis."
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options




#***************
# READ OPTIONS *
#***************

# read genes

if ( ! is.null(opt$genes) ) {
  
  if ( opt$genes == 'stdin' ) {
    
    G <- read.table( file = 'stdin', header = FALSE, quote = NULL)
    
  } else {
    
    G <- read.table( file = opt$genes, header = FALSE, quote = NULL)
    
  }
  
  colnames(G) <- "ids"
  
} else {
  
  print("Missing genes' list!")
  quit(save = 'no')
  
}


# read universe

if ( ! is.null(opt$universe) ) {
  
  if ( opt$universe == 'stdin' ) {
    
    U <- read.table( file = 'stdin', header = FALSE, quote = NULL)
    
  } else {
    
    U <- read.table( file = opt$universe, header = FALSE, quote = NULL)
    
  }
  
  colnames(U) <- "ids"
  
} else {
  
  print("Missing universe's list!")
  quit(save = 'no')
  
}


#************
# LIBRARIES *
#************

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("DBI"))
suppressPackageStartupMessages(library("GO.db"))
if (opt$species == "human") { suppressPackageStartupMessages(library("org.Hs.eg.db")) }
if (opt$species == "mouse") { suppressPackageStartupMessages(library("org.Mm.eg.db")) }
if (opt$species == "dmel") { suppressPackageStartupMessages(library("org.Dm.eg.db")) }
suppressPackageStartupMessages(library("GOstats"))
suppressPackageStartupMessages(library("plyr"))



#********
# BEGIN *
#********

if (opt$species == "human") {
  
  universe = unlist(mget(as.character(U$ids), org.Hs.egENSEMBL2EG, ifnotfound=NA))

} else if (opt$species == "mouse") {
  
  universe = unlist(mget(as.character(U$ids), org.Mm.egENSEMBL2EG, ifnotfound=NA))

} else if (opt$species == "dmel") {
  
  universe = unlist(mget(as.character(U$ids), org.Dm.egENSEMBL2EG, ifnotfound=NA))

}


createParams = function(x, species = "human") {
  
  if (species == "human") {
    
    ann = "org.Hs.eg.db"
    geneset = unlist(mget(x, org.Hs.egENSEMBL2EG, ifnotfound=NA))
  
  }
  
  if (species == "mouse") {
    
    ann = "org.Mm.eg.db"
    geneset = unlist(mget(x, org.Mm.egENSEMBL2EG, ifnotfound=NA))
  
  }
  
  if (species == "dmel") {
    
    ann = "org.Dm.eg.db"
    geneset = unlist(mget(x, org.Dm.egENSEMBL2EG, ifnotfound=NA))
  
  }
  
  sprintf("%s foreground genes; %s with a corresponding entrez id", length(x), length(unique(geneset)))
  
  pv = 1 
  params = new( "GOHyperGParams",
                geneIds = geneset,
                universeGeneIds = universe,
                annotation = ann,
                ontology = opt$category,
                pvalueCutoff = pv,
                conditional = TRUE,
                testDirection = 'over' )
  return(params) 

}

res = hyperGTest(createParams(as.character(G$ids), opt$species))

df = summary(res)
df$fdr = p.adjust(df$Pvalue, method = "BH")



# barplot of significant GO terms (fdr < 0.05)

sub.df <- df[df$fdr < opt$fdr_cutoff, ]

if (nrow(sub.df) > 0) {
  
  if (nrow(sub.df) <= 20) {
    
    p <- ggplot(sub.df,
                aes(x=reorder(Term, -log10(fdr + opt$pseudocount)), y=-log10(fdr + opt$pseudocount))) +
      geom_bar(stat="identity", fill="orange", colour="black") +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5, size=15),
            axis.text = element_text(size=opt$textsize),
            axis.title = element_text(size=15)) +
      coord_flip() +
      ylab("-log10(fdr)") +
      xlab(paste0("GO terms - ", opt$category))
    
    if ( ! is.null(opt$title) ) {
      
      p <- p + labs(title = opt$title)
      
    }
    
  } else {
    
    x <- 1:nrow(sub.df)
    sub.df$facet <- ceiling(x/20)
    my.n <- max(as.numeric(sub.df$facet))
    p <- list()
    
    for (i in 1:my.n) {
      
      p[[i]] <- ggplot(sub.df[sub.df$facet == i, ],
                         aes(x=reorder(Term, -log10(fdr + opt$pseudocount)), y=-log10(fdr + opt$pseudocount))) +
        geom_bar(stat="identity", fill="orange", colour="black") +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour = "black"),
              plot.title = element_text(hjust = 0.5, size=15),
              axis.text = element_text(size=opt$textsize),
              axis.title = element_text(size=15)) +
        coord_flip() +
        ylab("-log10(fdr)") +
        xlab(paste0("GO terms - ", opt$category))
      
      if ( ! is.null(opt$title) ) {
        
        p[[i]] <- p[[i]] + labs(title = opt$title)
        
      }
      
      
      
    }
    
    
    
  }
  

  
} else {
  
  p <- NULL
  
}

#*********
# OUTPUT *
#*********

output = sprintf("%s/%s.%s", opt$directory, opt$output, opt$category)
write.table(df, file=sprintf("%s.tsv", output), quote=F, sep="\t", row.names=F)
htmlReport(res, file=sprintf("%s.html", output))

if (! is.null(p)) {
  
  pdf(file = sprintf("%s.pdf", output), paper="a4r", width=32, height=14)
  if (nrow(sub.df) <= 20) {
    print(p)
  } else {
    for (i in 1:my.n) {
      print(p[i])
      
    }
    
  }
  dev.off()
}


#*******
# EXIT *
#*******

q( save = 'no' )

