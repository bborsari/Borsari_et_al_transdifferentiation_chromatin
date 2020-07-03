.libPaths("/nfs/users2/rg/bborsari/software/R-3.5.2/library")


#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option ( c("--output", "-o"), default = "out.6.groups.tsv",
                help = "Output file name. 'stdout' for standard output [default=%default]." ),
  
  make_option( c("--metadata", "-d"), 
               help = "Metadata .tsv file"),
  
  make_option( c("--expression_matrix", "-e"),
               help = "Expression matrix"),
  
  make_option( c("--mark_matrix", "-m"),
               help = "Mark matrix"),
  
  make_option( c("--alignment_matrix", "-a"),
               help = "Alignment matrix"),
  
  make_option( c("--coverage_matrix", "-C"),
               help = "Coverage matrix from Zerone"),
  
  make_option( c("--genes_intersecting_peaks", "-G"),
               help = "List of genes intersecting peaks"),
  
  make_option( c("--maSigPro_genes", "-s"),
               help = "List of genes significantly variable according to maSigPro")
  
)


parser <- OptionParser(
  usage = "%prog [options]", 
  option_list=option_list,
  description = "\nClassifies genes as positively correlated, negatively correlated, not correlated, stable, no_peak, peak_not_TSS."
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options


#*****************
# LOAD LIBRARIES *
#*****************

suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))


#********************
# DEBUGGING OPTIONS *
#********************

# setwd("H3K36me3/QN.merged")
# metadata <- read.table ( file = "metadata.tsv" , quote = NULL, header = T, sep="\t" )
# e.matrix <- read.table ( file = "expression.matrix.tsv", quote = NULL, header = T, sep="\t" )
# m.matrix <- read.table ( file = "H3K36me3.matrix.tsv", quote = NULL, header = T, sep="\t" )
# alignment.matrix <- read.table ( file = "NW.pr/H3K36me3.FDR.table.all.genes.tsv", quote = NULL, header = T, sep="\t",
#                                         stringsAsFactors = F)
# coverage.matrix <- read.table ( file = "all.tp.peak.length.QN.merged.tsv", quote = NULL, header = T, sep="\t" )
# G.matrix <- read.table ( file = "genes.intersecting.peaks.tsv", quote = NULL, header = F, sep="\t" )
# S.matrix <- read.table ( file = "H3K36me3.QN.merged.maSigPro.out.tsv", quote = NULL, header = T, sep="\t" )



#***************
# READ OPTIONS *
#***************


# 1. metadata file
metadata <- read.table ( file = opt$metadata, quote = NULL, header = T, sep="\t" )
stopifnot(nrow(metadata) == 8030)
metadata$gene_id <- rownames(metadata)

# 2. expression matrix
e.matrix <- read.table ( file = opt$expression_matrix, quote = NULL, header = T, sep="\t" )
e.matrix <- e.matrix[rownames(metadata), ]
stopifnot(identical(rownames(e.matrix), rownames(metadata)))


# 3. mark matrix
m.matrix <- read.table ( file = opt$mark_matrix, quote = NULL, header = T, sep="\t" )
m.matrix <- m.matrix[rownames(metadata), ]
stopifnot(identical(rownames(m.matrix), rownames(metadata)))


# 4. alignment matrix
alignment.matrix <- read.table ( file = opt$alignment_matrix, quote = NULL, header = T, sep="\t",
                                 stringsAsFactors = F )
alignment.matrix$gene_id <- rownames(alignment.matrix)


# 5. coverage matrix
coverage.matrix <- read.table ( file = opt$coverage_matrix, quote = NULL, header = T, sep="\t" )
coverage.matrix <- as.data.frame(as.matrix(coverage.matrix[, 1:12]) / coverage.matrix$gene_length)
coverage.matrix$sum_coverage <- apply(coverage.matrix, 1, sum)
coverage.matrix <- coverage.matrix[rownames(metadata), ]
stopifnot(identical(rownames(coverage.matrix), rownames(metadata)))
coverage.matrix$gene_id <- rownames(coverage.matrix)

# 6. genes intersecting peaks
G.matrix <- read.table ( file = opt$genes_intersecting_peaks, quote = NULL, header = F, sep="\t" )
colnames(G.matrix) <- "gene_id"
rownames(G.matrix) <- G.matrix$gene_id
G.matrix$peak <- "yes"
G.matrix <- G.matrix[rownames(G.matrix) %in% rownames(metadata), ]


# 7. genes significantly variable according to maSigPro
S.matrix <- read.table ( file = opt$maSigPro_genes, quote = NULL, header = T, sep="\t" )
S.matrix$gene_id <- rownames(S.matrix)
S.matrix <- S.matrix %>% separate(gene_id, c("gene", "id"), "\\.")
S.matrix$id <- NULL
colnames(S.matrix)[2] <- "gene_id"
rownames(S.matrix) <- S.matrix$gene_id
S.matrix <- S.matrix[rownames(S.matrix) %in% rownames(metadata), ]


# 8. output file name
output = ifelse(opt$output == "stdout", "", opt$output)




#********
# BEGIN *
#********

hours <- c(0,3,6,9,12,18,24,36,48,72,120,168)
colnames(m.matrix) <- hours
colnames(e.matrix) <- hours


###-----------------
### 1: compute cc
###-----------------

## Compute pearson cc between expression and mark level
df.mark.levels.pearson <- data.frame(gene_id = character(),
                                     correlation = numeric(),
                                     pvalue = numeric(),
                                     stringsAsFactors = F)

for (i in rownames(e.matrix)) {
  
  df.mark.levels.pearson[nrow(df.mark.levels.pearson)+1,] <- c(i, 
                                                               cor.test(t(e.matrix[i,]), t(m.matrix[i, ]))$estimate,
                                                               cor.test(t(e.matrix[i,]), t(m.matrix[i, ]))$p.value)
}

df.mark.levels.pearson$correlation <- as.numeric(df.mark.levels.pearson$correlation)
df.mark.levels.pearson$pvalue <- as.numeric(df.mark.levels.pearson$pvalue)
df.mark.levels.pearson$fdr <- p.adjust(df.mark.levels.pearson$pvalue, method = "fdr")


## Compute pearson cc between expression and mark coverage
df.mark.coverage.pearson <- data.frame(gene_id=character(),
                                       correlation = numeric(),
                                       pvalue = numeric(),
                                       stringsAsFactors = F)

for (i in rownames(e.matrix)) {
  
  df.mark.coverage.pearson[nrow(df.mark.coverage.pearson)+1,] <- c(i, 
                                                                   cor.test(t(e.matrix[i,]), t(coverage.matrix[i, 1:12]))$estimate,
                                                                   cor.test(t(e.matrix[i,]), t(coverage.matrix[i, 1:12]))$p.value)
}


df.mark.coverage.pearson$correlation <- as.numeric(df.mark.coverage.pearson$correlation)
df.mark.coverage.pearson$pvalue <- as.numeric(df.mark.coverage.pearson$pvalue)
df.mark.coverage.pearson$fdr <- p.adjust(df.mark.coverage.pearson$pvalue, method = "fdr")



###----------------------
### 2: merge dataframes
###----------------------


## add metadata class and presence/absence of peak in the region of interest
merged.df <- merge(x = metadata,
                   y = G.matrix, 
                   by = "gene_id",
                   all = T)
colnames(merged.df)[ncol(merged.df)] <- "presence_of_peak"
merged.df$presence_of_peak <- ifelse(is.na(merged.df$presence_of_peak), "no", merged.df$presence_of_peak)


## add sum of coverage over gene body across 12 tp
merged.df <- merge( x = merged.df,
                    y = coverage.matrix[, c("sum_coverage", "gene_id")],
                    by = "gene_id" )


## add whether the gene is significantly variable according to maSigPro
merged.df <- merge(x = merged.df,
                   y = S.matrix, 
                   by = "gene_id",
                   all = T)
colnames(merged.df)[ncol(merged.df)] <- "maSigPro"
merged.df$maSigPro <- ifelse(is.na(merged.df$maSigPro), "ns", "significant")


## add score, p-value and FDR of alignment
colnames(alignment.matrix) <- c("score_alignment", "pvalue_alignment", "fdr_alignment", "gene_id")
merged.df <- merge( x = merged.df,
                    y = alignment.matrix,
                    by = "gene_id",
                    all = T)


## add correlation, p-value and FDR of pearson cc between expression and mark level
colnames(df.mark.levels.pearson) <- c("gene_id",
                                      "pearson_estimate_mark_level",
                                      "pearson_pvalue_mark_level",
                                      "pearson_fdr_mark_level")
merged.df <- merge( x = merged.df,
                    y = df.mark.levels.pearson,
                    by = "gene_id" )


## add correlation, p-value and FDR of pearson cc between expression and mark coverage
colnames(df.mark.coverage.pearson) <- c("gene_id",
                                        "pearson_estimate_mark_coverage",
                                        "pearson_pvalue_mark_coverage",
                                        "pearson_fdr_mark_coverage")
merged.df <- merge( x = merged.df,
                    y = df.mark.coverage.pearson,
                    by = "gene_id" )

rownames(merged.df) <- merged.df$gene_id


###--------------------------
### 3: start classification
###--------------------------


## genes w/o a peak in the region of interest at any point of the process
no.peak.genes <- merged.df[(merged.df$presence_of_peak == "no" & merged.df$sum_coverage == 0), "gene_id"]
if (length(no.peak.genes) > 0) {
  
  no.peak.df <- data.frame(gene_id = no.peak.genes,
                           group = "no_peak")
  
} else {
  
  no.peak.df <- data.frame(gene_id = c(),
                           group = c()) 
  
}

## genes w/o a peak in the region of interest but positive coverage over the gene body
peak.not.TSS.genes <- merged.df[(merged.df$presence_of_peak == "no" & merged.df$sum_coverage > 0), "gene_id"]
if (length(peak.not.TSS.genes) > 0) {
  
  peak.not.TSS.df <- data.frame(gene_id = peak.not.TSS.genes,
                                group = "peak_not_TSS")
  
} else {
  
  peak.not.TSS.df <- data.frame(gene_id = c(),
                                group = c()) 
  
}



## genes w/ a peak in the region of interest but absence of significant changes according to maSigPro
stable.genes <- merged.df[(merged.df$presence_of_peak == "yes" & merged.df$maSigPro == "ns"), "gene_id"]
stable.df <- data.frame(gene_id = stable.genes,
                        group = "stable")

## genes left to classify
tmp.group <- c()
tmp <- merged.df[merged.df$presence_of_peak == "yes" & merged.df$maSigPro == "significant", ]
for ( r in rownames(tmp) ) {
  
  if ( ( !is.na(tmp[r, "pearson_fdr_mark_level"]) & tmp[r, "pearson_fdr_mark_level"] < 0.05) & # cc with mark level significant
       ( !is.na(tmp[r, "pearson_fdr_mark_coverage"]) & tmp[r, "pearson_fdr_mark_coverage"] < 0.05) & # cc with mark coverage significant
       ( !is.na(tmp[r, "fdr_alignment"]) & tmp[r, "fdr_alignment"] < 0.05) & # significant alignment
       ( !is.na(tmp[r, "pearson_estimate_mark_level"]) & tmp[r, "pearson_estimate_mark_level"] >= 0.6) & # cc with mark level positive
       ( !is.na(tmp[r, "pearson_estimate_mark_coverage"]) & tmp[r, "pearson_estimate_mark_coverage"] >= 0.6) # cc with mark coverage positive
  ) {
    
    group = "positively_correlated"
    
  } else if ( (!is.na(tmp[r, "pearson_fdr_mark_level"]) & tmp[r, "pearson_fdr_mark_level"] < 0.05) & # cc with mark level significant
              (!is.na(tmp[r, "pearson_fdr_mark_coverage"]) & tmp[r, "pearson_fdr_mark_coverage"] < 0.05) & # cc with mark coverage significant
              (!is.na(tmp[r, "fdr_alignment"]) & tmp[r, "fdr_alignment"] >= 0.05) & # not significant alignment
              (!is.na(tmp[r, "pearson_estimate_mark_level"]) & tmp[r, "pearson_estimate_mark_level"] <= -0.6) & # cc with mark level negative
              (!is.na(tmp[r, "pearson_estimate_mark_coverage"]) & tmp[r, "pearson_estimate_mark_coverage"] <= -0.6) # cc with mark coverage negative 
  ) {
    
    group = "negatively_correlated"
    
  } else if ( (!is.na(tmp[r, "pearson_fdr_mark_level"]) & tmp[r, "pearson_fdr_mark_level"] < 0.05) & # cc with mark level significant 
              (!is.na(tmp[r, "pearson_fdr_mark_coverage"]) & tmp[r, "pearson_fdr_mark_coverage"] < 0.05) & # cc with mark coverage significant
              (!is.na(tmp[r, "pearson_estimate_mark_level"]) & tmp[r, "pearson_estimate_mark_level"] >= 0.6) & # cc with mark level positive
              (!is.na(tmp[r, "pearson_estimate_mark_coverage"]) & tmp[r, "pearson_estimate_mark_coverage"] >= 0.6) # cc with mark coverage positive
  ) {
    
    group = "positively_correlated"
    
  } else if ( (!is.na(tmp[r, "pearson_fdr_mark_level"]) & tmp[r, "pearson_fdr_mark_level"] < 0.05) &  # cc with mark level significant
              (!is.na(tmp[r, "fdr_alignment"]) & tmp[r, "fdr_alignment"] < 0.05) & # alignment significant
              (!is.na(tmp[r, "pearson_estimate_mark_level"]) & tmp[r, "pearson_estimate_mark_level"] >= 0.6) ) {
    
    group="positively_correlated"
    
  } else if ( (!is.na(tmp[r, "pearson_fdr_mark_coverage"]) & tmp[r, "pearson_fdr_mark_coverage"] < 0.05) & # cc with mark coverage significant
              (!is.na(tmp[r, "fdr_alignment"]) & tmp[r, "fdr_alignment"] < 0.05) & # alignment significant
              (!is.na(tmp[r, "pearson_estimate_mark_coverage"]) & tmp[r, "pearson_estimate_mark_coverage"] >= 0.6) # cc with mark coverage positive
  ) {
    
    group="positively_correlated"
    
  } else if ( (!is.na(tmp[r, "pearson_fdr_mark_level"]) & tmp[r, "pearson_fdr_mark_level"] < 0.05) & # cc with mark level significant
              (!is.na(tmp[r, "pearson_fdr_mark_coverage"]) & tmp[r, "pearson_fdr_mark_coverage"] < 0.05) & # cc with mark coverage significant
              (!is.na(tmp[r, "pearson_estimate_mark_level"]) & tmp[r, "pearson_estimate_mark_level"] <= -0.6) & # cc with mark level negative
              (!is.na(tmp[r, "pearson_estimate_mark_coverage"]) & tmp[r, "pearson_estimate_mark_coverage"] <= -0.6) # cc with mark coverage negative
  ){
    
    group="negatively_correlated"
    
  } else if ( (!is.na(tmp[r, "pearson_fdr_mark_coverage"]) & tmp[r, "pearson_fdr_mark_coverage"] < 0.05) & # cc coverage significant & negative 
              (!is.na(tmp[r, "fdr_alignment"]) & tmp[r, "fdr_alignment"] >= 0.05) & # alignment not significant 
              (!is.na(tmp[r, "pearson_estimate_mark_coverage"]) & tmp[r, "pearson_estimate_mark_coverage"] <= -0.6) # cc with mark coverage negative
  ){
    
    group="negatively_correlated"
    
  } else if ( (!is.na(tmp[r, "pearson_fdr_mark_level"]) & tmp[r, "pearson_fdr_mark_level"] < 0.05) &  # cc with mark level significant
              (!is.na(tmp[r, "fdr_alignment"]) & tmp[r, "fdr_alignment"] >= 0.05) & # alignment not significant
              (!is.na(tmp[r, "pearson_estimate_mark_level"]) & tmp[r, "pearson_estimate_mark_level"] <= -0.6) # cc with mark level negative
  ) {
    
    group="negatively_correlated"
    
  } else {group = "not_correlated"}
  
  tmp.group <- c(tmp.group, group)
  
}

variable.genes <- rownames(tmp)
variable.df <- data.frame(gene_id = variable.genes,
                          group = tmp.group)

all.genes.groups <- rbind(no.peak.df,
                          peak.not.TSS.df,
                          stable.df,
                          variable.df)

merged.df <- merge( x = merged.df,
                    y = all.genes.groups,
                    by = "gene_id" )

rownames(merged.df) <- merged.df$gene_id
merged.df$gene_id <- NULL
merged.df <- merged.df[rownames(metadata), ]
stopifnot(identical(rownames(merged.df), rownames(metadata)))

#*********
# OUTPUT *
#*********

write.table( merged.df, output, quote = FALSE, sep='\t', row.names = T )


#*******
# EXIT *
#*******


quit( save = "no" )



