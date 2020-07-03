.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option ( c("--input", "-i"),
                help = "matrix of expression and chromatin profiles" ),
  
  make_option ( c("--header"), default = TRUE,
                help = "Whether the input matrix has a header [default=%default]." ),
  
  make_option ( c("--degree", "-d"), type = 'numeric', default = 2,
                help = "The regression degree [default=%default]." ),
  
  make_option ( c("--output", "-o"), default = "out.loess.pr.tsv",
                help = "Output file name. 'stdout' for standard output [default=%default]." )
  
)


parser <- OptionParser(
  usage = "%prog [options]", 
  option_list=option_list,
  description = "\nPerforms polynomial regression on each row of the input matrix."
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options




#***************
# READ OPTIONS *
#***************


if ( !is.null(opt$input) ) {
  
  input.matrix <- read.table( file = opt$input, header = opt$header, quote = NULL, sep="\t" )

} else {
  
  print('Missing input!')
  quit(save = 'no')
  
}


output = ifelse(opt$output == "stdout", "", opt$output)


if ( opt$header == FALSE ) {
  
  rownames(input.matrix) <- input.matrix$V1
  input.matrix$V1 <- NULL
  
}



#********
# BEGIN *
#********

input.matrix.loess <- data.frame(stringsAsFactors = F)


for (r in 1:nrow(input.matrix)) {

  # loess on row r
  df <- as.data.frame(t(input.matrix[r,]))
  colnames(df) <- "level"
  df$time <- c(0,3,6,9,12,18,24,36,48,72,120,168)
  loess.vector <- predict(loess(level~time, data = df, degree = opt$degree), 
                          0:168)[c(0,3,6,9,12,18,24,36,48,72,120,168)+1]

  input.matrix.loess <- rbind(input.matrix.loess,
                              loess.vector)

}

rownames(input.matrix.loess) <- rownames(input.matrix)
colnames(input.matrix.loess) <- colnames(input.matrix)




#*********
# OUTPUT *
#*********


write.table(input.matrix.loess, file=output, sep="\t", quote=F)



#*******
# EXIT *
#*******

quit(save="no")
