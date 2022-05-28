.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")



#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  
  make_option( c ("-i", "--input"), default = 'stdin',
               help = "input bed file [default = %default]." ),
  
  
  make_option( c("-o", "--output"), default = 'out.tsv', 
               help = "output file")
  
)


parser <- OptionParser(
  usage = "%prog [options] files", 
  option_list=option_list,
  description = "\nWrapper for HMM"
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options


#********
# BEGIN *
#********

# 1. read bed file
if ( ! is.null(opt$input) ) {
  
  if ( opt$input == 'stdin' ) {
    
    m <- read.table( file = 'stdin', header = FALSE, quote = NULL, sep="\t")
    
  } else {
    
    m <- read.table( file = opt$input, header = FALSE, quote = NULL, sep="\t")
    
  }
  
  colnames(m) <- c("chrom", "start", "end", "value")
  
} else {
  
  print("Missing input file!")
  quit(save = 'no')
  
}


# 2. read chrom sizes
chrom_sizes <- read.table("~/projects/ERC/human/ERC_human_hub/chrom.sizes.tsv",
                          h=F, sep="\t")
colnames(chrom_sizes) <- c("chrom", "length")


# 3. fix wrong peak coordinates
m <- merge(m, chrom_sizes, by = "chrom")
m$end <- ifelse(m$end <= m$length, m$end, m$length)
m <- m[, 1:4]


# 4. save file
output = opt$output
write.table(m, file=opt$output, quote=F, sep="\t", row.names=F, col.names = F)
quit(save="no")
