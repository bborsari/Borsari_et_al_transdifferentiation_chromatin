.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option ( c("--input_matrix", "-i"), 
                help = "Input z-score matrix (after polynomial regression)." ),
  make_option ( c("--output_name", "-o"),
                help = "Output filename." )
  
)

parser <- OptionParser(
  usage = "%prog [options] files", 
  option_list=option_list)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options


#************
# FUNCTIONS *
#************

rescale <- function(x){
  return((x -min(x))/(max(x) - min(x)))
}




#********
# BEGIN *
#********


# 1. read input matrix
if (!( is.null(opt$input_matrix) )) {
  
  m <- read.table(opt$input_matrix, h=T, sep="\t")
  
} else {
  
  print("Missing input matrix!")
  quit(save="no")
  
}


# 2. rescale matrix rows to range 0-1 
m.rescaled <- t(apply(m, 1, rescale))


# 3. save output
write.table(m.rescaled, file = opt$output_name, col.names = T, row.names = T, 
            quote = F, sep="\t")

