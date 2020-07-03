#!/usr/bin/env Rscript



#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option( c("--input", "-i"), default = 'stdin',
               help = "Input matrix. 'stdin' for reading from standard input [default = %default]." ),
  
  make_option( c("--output", "-o"), default = 'out.p-adjust.tsv',
               help = "Output file name. 'stdout' for standard output [default = %default]." ),
  
  make_option( c("--header"), default = TRUE,
               help = "Whether the input matrix has a header [default = %default]." ),
  
  make_option( c("--column", "-c"), type = 'numeric', default = 1,
               help = "The column that corresponds to the p-value [default = %default]." ),
  
  make_option( c("--method", "-m"), default = 'fdr',
               help = "holm, hochberg, hommel, bonferroni, BH, BY, fdr [default = %default]." ),
  
  make_option( c("--round", "-r"), default = NULL, type = 'numeric',
               help = "The number of decimal digits to keep for the adjusted p-values [default = %default].")
  
)


parser <- OptionParser(
  usage = "%prog [options] file", 
  option_list=option_list,
  description = "\nProvided an input matrix, it attaches an additional column with the corrected p-value."
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options





#***************
# READ OPTIONS *
#***************


if ( ! is.null(opt$input) ) {
  
  if ( opt$input == 'stdin' ) {
    
    input.matrix <- read.table( file = 'stdin', header = opt$header, quote = NULL )
    
  } else {
    
    input.matrix <- read.table( file = opt$input, header = opt$header, quote = NULL )
    
  }
  
} else {
  
  print('Missing input matrix!')
  quit(save = 'no')
  
}


output = ifelse(opt$output == "stdout", "", opt$output)



#********
# BEGIN *
#********


if ( opt$header == FALSE ) {
  
  rownames(input.matrix) <- input.matrix$V1
  input.matrix$V1 <- NULL
  
}


input.matrix$adjusted <- p.adjust(input.matrix[, opt$column], method = opt$method)

if (!is.null(opt$round)) {
  
  input.matrix$adjusted <- round(input.matrix$adjusted, opt$round)
  
}

colnames(input.matrix)[ncol(input.matrix)] <- opt$method

write.table( input.matrix, output, quote = FALSE, sep='\t' )




#*******
# EXIT *
#*******


quit( save = "no" )

