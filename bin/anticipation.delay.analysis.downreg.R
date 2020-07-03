.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option ( c("--input_matrix", "-i"), 
                help = "Input matrix (after polynomial regression)." ),
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


my.function <- function(my.df.rescaled) {
  
  my.res.df <- data.frame(stringsAsFactors = F)
  
  for ( r in rownames(my.df.rescaled) ) {
    
    current.gene.df <- as.data.frame(my.df.rescaled[r, ])
    colnames(current.gene.df) <- "values"
    current.gene.df$perc_75 <- abs(0.75 - current.gene.df$values)
    current.gene.df$perc_50 <- abs(0.5 - current.gene.df$values)
    current.gene.df$perc_25 <- abs(0.25 - current.gene.df$values)
    current.gene.df$perc_0 <- abs(0 - current.gene.df$values)
    
    current.gene.df$tp <- 1:12
    
    tmp.df <- data.frame(perc_75 = current.gene.df[order(current.gene.df$perc_75), "tp"],
                         perc_50 = current.gene.df[order(current.gene.df$perc_50), "tp"],
                         perc_25 = current.gene.df[order(current.gene.df$perc_25), "tp"],
                         perc_0 = current.gene.df[order(current.gene.df$perc_0), "tp"])
    my.v <- c()
    
    my.keep = T
    
    while (my.keep) {
      
      # 75 vs. 50%
      for (i in 1:nrow(tmp.df)) {
        
        if ( (as.numeric(tmp.df[i, 1]) <= as.numeric(tmp.df[1, 2])) &
             (as.numeric(tmp.df[i, 1]) <= as.numeric(tmp.df[1, 3])) & 
             (as.numeric(tmp.df[i, 1]) <= as.numeric(tmp.df[1, 4])) )  {
          
          my.v <- c(my.v, tmp.df[i, 1])
          break
          
        }
        
      }
      
      
      # 50 vs. 25%
      for (i in 1:nrow(tmp.df)) {
        
        if ( as.numeric(tmp.df[i, 2]) <= as.numeric(tmp.df[1, 3]) &
             (as.numeric(tmp.df[i, 2]) <= as.numeric(tmp.df[1, 4])) ) {
          
          my.v <- c(my.v, tmp.df[i, 2])
          break
          
        }
        
      }
      
      
      # 25 vs. 0%
      for (i in 1:nrow(tmp.df)) {
        
        if (as.numeric(tmp.df[i, 3]) <= as.numeric(tmp.df[1, 4])) {
          
          my.v <- c(my.v, tmp.df[i, 3], tmp.df[1, 4])
          break
          
        }
        
      }
      
      names(my.v) <- 1:4
      
      if (sum(names(sort(my.v)) == 1:4) == 4) {
        
        my.keep = F
        
      } else {
        
        tmp.df[1, ] <- my.v
        my.v <- c()
        
      }
      
      
    }
    
    
    my.res.df <- rbind(my.res.df, my.v)
    
  }
  
  return(my.res.df)
  
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


# 3. get matrix of changes for 75, 50, 25 and 0% of downregulation
m.out <- my.function(my.df.rescaled = m.rescaled)
colnames(m.out) <- c("perc_75", "perc_50", "perc_25", "perc_0")
rownames(m.out) <- rownames(m)

# 4. save output
write.table(m.out, file = opt$output_name, col.names = T, row.names = T, 
            quote = F, sep="\t")



