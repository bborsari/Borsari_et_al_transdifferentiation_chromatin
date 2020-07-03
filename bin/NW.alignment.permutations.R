.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option ( c("--input", "-i"),
                help = "matrix of expression and chromatin profiles" ),
  
  make_option ( c("--batch", "-b"), type = "numeric",
                help = "batch number" ),
  
  make_option ( c("--dir", "-d"), default = ".",
                help = "output directory [default=%default]" ),
  
  make_option ( c("--header"), default = TRUE,
                help = "Whether the input matrix has a header [default=%default]." ),
  
  make_option ( c("--t1"), type = "numeric",
                help = "number of time points profile #1" ),
  
  make_option ( c("--t2"), type = "numeric",
                help = "number of time points profile #2" ),
  
  make_option ( c("--replicate", "-r"), default = NULL,
                help = "sample replicate [default=%default]" ),
  
  make_option ( c("--permutations", "-p"), default = 250000, type = "numeric",
                help = "number of permutations [default=%default]" )

)


parser <- OptionParser(
  usage = "%prog [options]", 
  option_list=option_list,
  description = "\nProvided a series of expression and chromatin profiles, for each pair performs a NW alignment on the original profiles and a number of alignments on permutations of the two profiles (default = 250000). It returns the distance of the alignments (original and permuted)."
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options


#***************
# READ OPTIONS *
#***************


if ( !is.null(opt$batch) ) {
  
  batch <- opt$batch

} else {
  
  print( "Bad inline parameters" )
  quit( save = "no" )
  
}

if ( !is.null(opt$input) ) {
  
  data_matrix <- read.table ( file = opt$input, quote = NULL, header = opt$header )

} else {
  
  print( "Missing input!" )
  quit( save = "no" )
  
}


if ( !is.null(opt$t1) ) {

  t1 <- opt$t1

} else {
  
  print( 'Missing t1!' )
  quit( save = "no" )
  
}


if ( !is.null( opt$t2) ) {
  
  t2 <- opt$t2
  
} else {
  
  print( 'Missing t2!' )
  quit( save = "no" )
  
}






#*******************
# DEFINE FUNCTIONS *
#*******************


get_distances <- function ( expression_sequence, mark_sequence ) {
  
  dist_matrix <- matrix( nrow = (t1+1), ncol = (t2+1) )
  dist_matrix[1,1] <- 0
  
  for ( mark in 2:(t2+1) ) {
    
    for ( expression in 2:(t1+1) ) {
      
      dist_matrix[ expression,mark ] <- 
        abs( expression_sequence[expression] - mark_sequence[mark] )
      
    }
  
  }
  
  for ( mark in 2:(t2+1) ) { 
    
    dist_matrix[ 1,mark ] <- dist_matrix[ 2,mark ] 
    
  }
  
  for ( expression in 2:(t1+1) ) { 
    
    dist_matrix[ expression,1 ] <- dist_matrix[ expression,2 ] 
    
  }
  
  return(dist_matrix)
  
}


align_sequences <- function ( expression_sequence, mark_sequence ) {
  
  dist_matrix <- get_distances ( expression_sequence, mark_sequence )
  
  align_matrix <- matrix( nrow = (t1+1), ncol = (t2+1) )
  
  align_matrix[1,1] <- 0
  
  for ( mark in 2:(t2+1) ) {
    
    align_matrix[ 1,mark ]  <- dist_matrix[ 2,mark ] + align_matrix[ 1,(mark-1) ]
  
  }
  
  for ( expression in 2:(t1+1) ) {
    
    align_matrix[ expression,1 ] <- 
      dist_matrix[ expression,2 ] + align_matrix[ (expression-1),1 ]
  }
  
  for ( mark in 2:(t2+1) ) {
    
    for( expression in 2:(t1+1) ) {
      
      diag <- align_matrix[ (expression-1),(mark-1) ]
      up <- align_matrix[ (expression-1),mark ]
      left <- align_matrix[ expression,(mark-1) ]
      
      align_matrix[ expression,mark ] <- 
        min( c( diag, up, left) ) + dist_matrix[ expression,mark ]
    
    }
  
  }

  return( align_matrix[ (t1+1),(t2+1) ] )

}


#********
# BEGIN *
#********



batch_size <- 50

pair_start <- (batch-1)*batch_size + 1
pair_end <- batch*batch_size

if ( opt$header == FALSE ) {
  
  rownames(data_matrix) <- data_matrix$V1
  data_matrix$V1 <- NULL
  colnames(data_matrix) <- c( 'expression_000h',
                              'expression_003h',
                              'expression_006h',
                              'expression_009h',
                              'expression_012h',
                              'expression_018h',
                              'expression_024h',
                              'expression_036h',
                              'expression_048h',
                              'expression_072h',
                              'expression_120h',
                              'expression_168h',
                              'mark_000h',
                              'mark_003h',
                              'mark_006h',
                              'mark_009h',
                              'mark_012h',
                              'mark_018h',
                              'mark_024h',
                              'mark_036h',
                              'mark_048h',
                              'mark_072h',
                              'mark_120h',
                              'mark_168h' )
  
}


if ( !is.null( opt$replicate) ) {
  
  replicate <- opt$replicate
  file_name <- paste0( opt$dir, "/", 'Batch.', as.character(batch), '.data.', replicate )
  
} else {
  
  file_name <- paste0( opt$dir, "/", 'Batch.', as.character(batch), '.data' )
  
}


write( paste( c( 'score','permutations' ), collapse='\t' ), file=file_name )

options( width=150 )

data_matrix$base <- 0

data_matrix <- as.matrix( data_matrix )

all_names <- rownames( data_matrix )


expression_cols <- c( 'base',
                      'expression_000h',
                      'expression_003h',
                      'expression_006h',
                      'expression_009h',
                      'expression_012h',
                      'expression_018h',
                      'expression_024h',
                      'expression_036h',
                      'expression_048h',
                      'expression_072h',
                      'expression_120h',
                      'expression_168h' )

mark_cols <- c( 'base',
                'mark_000h',
                'mark_003h',
                'mark_006h',
                'mark_009h',
                'mark_012h',
                'mark_018h',
                'mark_024h',
                'mark_036h',
                'mark_048h',
                'mark_072h',
                'mark_120h',
                'mark_168h' )



for ( index in pair_start:pair_end ) {
  
  if ( index > nrow(data_matrix) ) {
   
    break
     
  }

  gene_index <- all_names[index]
  
  expression_sequence <- data_matrix[ gene_index,expression_cols ]
  mark_sequence <- data_matrix[ gene_index,mark_cols ]
  current_score  <- sprintf("%.6f", align_sequences (expression_sequence,mark_sequence) )
  
  observations <- vector( length = opt$permutations )
  
  print(gene_index)
  
  for ( perm in 1:opt$permutations ) {
    
    set.seed(perm)
    expression_perm <- c("base" = 0)
    expression_perm <- c(expression_perm, sample( expression_sequence[2:length(expression_sequence)] ))
    
    mark_perm <- c("base" = 0)
    mark_perm <- c(mark_perm, sample( mark_sequence[2:length(mark_sequence)] ))

    print(perm)
    print(expression_perm)
    print(mark_perm)
    
    observations[perm] <- sprintf("%.6f", align_sequences (expression_perm,mark_perm) )
  
  }
  
  observations <- paste( sort( observations ), collapse = ',' )
  
  write( paste( c( gene_index,
                   current_score,
                   observations), 
                collapse='\t'),
         file=file_name, append = T )

}



#*******
# EXIT *
#*******


quit( save = "no" )
