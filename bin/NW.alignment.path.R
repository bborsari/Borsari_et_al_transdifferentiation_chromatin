#!/usr/bin/env Rscript



#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option ( c("--input", "-i"),
                help = "matrix of expression and chromatin profiles" ),
  
  make_option ( c("--output", "-o"), default = "out.NW.path.tsv",
                help = "Output file name. 'stdout' for standard output [default=%default]." ),
  
  make_option ( c("--header"), default = TRUE,
                help = "Whether the input matrix has a header [default=%default]." ),
  
  make_option ( c("--t1"), type = "numeric",
                help = "number of time points profile #1" ),
  
  make_option ( c("--t2"), type = "numeric",
                help="number of time points profile #2" )
  
)


parser <- OptionParser(
  usage = "%prog [options]", 
  option_list=option_list,
  description = "\nProvided a series of expression and chromatin profiles, for each pair performs a NW alignment on the original profiles and reconstructs the path. For each step in the alignment reports either 'match' (diagonal movement) or 'insertion' (up/left movement). Please consider to run NW.bidirectional.matches.py on the output of this script."
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options



#***************
# READ OPTIONS *
#***************


if ( !is.null(opt$input) ) {
  
  data_matrix <- read.table ( file = opt$input, quote = NULL, header = opt$header )
  
} else {
  
  print('Missing input!')
  quit(save = 'no')
  
}


if ( !is.null(opt$t1) ) {
  
  t1 <- opt$t1
  
} else {
  
  print('Missing t1!')
  quit(save = 'no')
  
}


if ( !is.null(opt$t2) ) {
  
  t2 <- opt$t2
  
} else {
  
  print('Missing t2!')
  quit(save = 'no')
  
}

output = ifelse(opt$output == "stdout", "", opt$output)


#*******************
# DEFINE FUNCTIONS *
#*******************


get_distances <- function( expression_sequence, mark_sequence ) {
  
  dist_matrix  <- matrix( nrow = (t1+1), ncol = (t2+1) )
  dist_matrix[1,1] <- 0
  
  for ( mark in 2:(t2+1) ) {
    
    for( expression in 2:(t1+1) ) {
      
      dist_matrix[ expression,mark ] <- 
        abs( expression_sequence [ expression ] - mark_sequence [ mark ] )
    
    }
    
  }
  
  for ( mark  in 2:(t2+1) ) { 
  
      dist_matrix[ 1,mark ] <- dist_matrix[ 2,mark ] 
  }
  
  for ( expression in 2:(t1+1) ) { 
    
    dist_matrix[ expression,1 ] <- dist_matrix[ expression,2 ] 

  }
  
  return(dist_matrix)
  
}




align_sequences <- function ( expression_sequence,
                              expression_cols,
                              mark_sequence,
                              mark_cols,
                              gene_index,
                              output ) {

    dist_matrix <- get_distances ( expression_sequence, mark_sequence )

    align_matrix <- matrix ( nrow = (t1+1), ncol = (t2+1) )
    
    align_matrix [ 1,1 ] <- 0
    
    for ( mark in 2:(t2+1) ) {
      
      align_matrix[ 1,mark ] <- dist_matrix [ 2,mark ] + align_matrix [ 1,(mark-1) ]
      # we add dist_matrix[ expression, 2 ] because we don't allow initial deletions
    
    }

    for ( expression in 2:(t1+1) ) {
      
      align_matrix[ expression,1 ] <- 
        dist_matrix[ expression,2 ] + align_matrix[ (expression-1),1 ] 
      # we add dist_matrix[ expression, 2 ] because we don't allow initial deletions
    
    }

    for ( mark in 2:(t2+1) ) {
      
      for ( expression in 2:(t1+1) ) {
        
        diag <- align_matrix[ (expression-1),(mark-1) ]
        up <- align_matrix[ (expression-1),mark ]
        left <- align_matrix[ expression,(mark-1) ]
        
	    align_matrix[ expression,mark ] <- 
	      min(c(diag, up, left)) + dist_matrix[ expression, mark ]
	    # min distance accumulated until this point from whatever direction + 
	    # distance of the current pair
	    
	    }
    
    }
    
    expression <- (t1+1)
    mark <- (t2+1)

    while ( expression > 1 && mark > 1 ) {
      
      diag <- align_matrix[ (expression-1),(mark-1) ]
      up <- align_matrix[ (expression-1),mark ]
      left <- align_matrix[ expression,(mark-1) ]
      
      if( diag <= up && diag <= left ) {
      # giving priority to 'diag' steps is a way to penalize gaps (up or left) 
        
        pair <- paste( c( gene_index,
                          'matching',
                          diag,
                          expression_cols[expression],
                          mark_cols[mark],
                          expression_sequence[expression],
                          mark_sequence[mark],
                          dist_matrix[ expression,mark ]),
                       collapse = "\t")
        
        write(pair, output, append = T)
        expression <- expression - 1
        mark <- mark - 1
        next
        
      }
      
      if ( up == left ) {
        
        write(paste(c(gene_index, expression, mark, up, left), collapse="\t"), 
              file="up.left.tsv", append = T)
        
      }
      
      if ( up <= diag && up <= left ) {
        
        pair <- paste( c( gene_index,
                          'insertion',
                          up,
                          expression_cols[expression],
                          mark_cols[mark],
                          expression_sequence[expression],
                          mark_sequence[mark],
                          dist_matrix[ expression,mark ] ), 
                       collapse = "\t")
        write(pair, output, append = T)
        expression <- expression - 1
        next
      
      }
      
      if ( left <= diag && left <= up ) {
        
        pair <- paste( c( gene_index,
                          'insertion',
                          left,
                          expression_cols[expression],
                          mark_cols[mark],
                          expression_sequence[expression],
                          mark_sequence[mark],
                          dist_matrix[ expression,mark ] ),
                       collapse = "\t")
        write(pair, output, append = T)
        mark <- mark - 1
        next
      }
      
    }
    
    while( expression > 1 ) {
      
      up  <- align_matrix[ (expression-1),1 ]
      pair  <- paste( c ( gene_index,
                          'terminal-expression',
                          up,
                          expression_cols[expression],
                          mark_cols[2],
                          expression_sequence[expression],
                          mark_sequence[2],
                          dist_matrix[ expression,2 ] ),
                      collapse = "\t" )
      write(pair, output, append = T)
      expression <- expression - 1

    }

    while ( mark > 1 ) {
      
      left  <- align_matrix[ 1,(mark-1) ]
      pair  <- paste( c ( gene_index,
                          'terminal-mark',
                          left,
                          expression_cols[2],
                          mark_cols[mark],
                          expression_sequence[2],
                          mark_sequence[mark],
                          dist_matrix[ 2,mark ] ),
                      collapse = "\t")
      write(pair, output, append = T)
      mark <- mark - 1
    
    }

}
            




#********
# BEGIN *
#********


options(width=150)


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


data_matrix$base <- 0
data_matrix <- as.matrix(data_matrix)
all_names <- rownames(data_matrix)

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

write ( paste (c( 'id',
                  'step-type',
                  'score',
                  'expression_time_point',
                  'mark_time_point',
                  'expression_z',
                  'mark_z',
                  'distance'), 
               collapse='\t'), 
        file=output )



for ( index in 1:nrow (data_matrix) ) {
  
  gene_index <- all_names[index]
  
  expression_sequence <- data_matrix[ gene_index,expression_cols ]
  mark_sequence <- data_matrix[ gene_index,mark_cols ]
  align_sequences ( expression_sequence, 
                    expression_cols, 
                    mark_sequence,
                    mark_cols,
                    gene_index,
                    output )

}



#*******
# EXIT *
#*******


quit( save = "no" )

