.libPaths("/nfs/users2/rg/bborsari/software/R-3.6.3/library")



#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option( c("-s", "--states"), type = "numeric",
               help = "Number of states." ),
  
  make_option( c("-o", "--outFolder"), 
               help = "Output folder."),
  
  make_option( c("-e", "--expression"), default = T,
               help = "Whether to include expression data in the model.")
  
  
)


parser <- OptionParser(
  usage = "%prog [options] files", 
  option_list=option_list,
  description = "\nWrapper for HMM"
  
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options




#************
# LIBRARIES *
#************

library(reshape2)
library(depmixS4)



#********
# BEGIN *
#********


# 1. the marks we're analyzing
marks <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3",
           "H3K36me3", "H4K20me1", "H3K9me3", "H3K27me3")


# 2. prepare input matrix for multivariate HMM

## 2.1. read expression matrix
m <- read.table("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/expression/QN.merged/selected.genes.rep.2.3.after.QN.merged.tsv",
                h=T, sep="\t")


# 2.2 add random noise to
# expression levels of not expressed genes
for ( i in 10697:12248 ) {
  
  set.seed(i)
  m[i, ] <- rnorm(12, mean = 0, sd = 10^-16)
  
}


m$gene_id <- rownames(m)
genes <- m$gene_id
m <- melt(m)
m$gene_id <- paste(m$gene_id, m$variable, sep="_")
m$variable <- NULL
colnames(m)[ncol(m)] <- "expression"


## 2.3. read marks' matrices
for ( i in 1:9 ) {
  
  tmp <- read.table(paste0("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/",
                           marks[i], "/QN.merged/",
                           marks[i], ".matrix.after.QN.merged.tsv"), h=T, sep="\t")
  
  tmp$gene_id <- rownames(tmp)
  tmp <- melt(tmp)
  tmp$gene_id <- paste(tmp$gene_id, tmp$variable, sep="_")
  tmp$variable <- NULL
  colnames(tmp)[ncol(tmp)] <- marks[i]
  
  m <- merge(m, tmp, by = "gene_id")
  
  
}


rownames(m) <- m$gene_id
m$gene_id <- NULL



# 3. run HMM
set.seed(123)

if (opt$expression) {
  
  mod1 <- depmix(list(expression~1, H3K4me1~1, H3K4me2~1, H3K27ac~1, 
                      H3K9ac~1, H3K4me3~1, H3K36me3~1, 
                      H4K20me1~1, H3K9me3~1, H3K27me3~1), 
                 data = m, nstates = opt$states,
                 family = list(gaussian(), gaussian(), gaussian(), gaussian(), 
                               gaussian(), gaussian(), gaussian(),
                               gaussian(), gaussian(), gaussian()),
                 ntimes = rep(12, 12248))

} else {
  
  mod1 <- depmix(list(H3K4me1~1, H3K4me2~1, H3K27ac~1, 
                      H3K9ac~1, H3K4me3~1, H3K36me3~1, 
                      H4K20me1~1, H3K9me3~1, H3K27me3~1), 
                 data = m[, 2:10], nstates = opt$states,
                 family = list(gaussian(), gaussian(), gaussian(), 
                               gaussian(), gaussian(), gaussian(),
                               gaussian(), gaussian(), gaussian()),
                 ntimes = rep(12, 12248))

  
}


fm1 <- fit(mod1, verbose = FALSE, emc=em.control(rand=T))



# 4. save df with states sequences
fm1.df <- posterior(fm1)
x <- matrix(fm1.df$state, byrow=T, nrow=12248)
x <- as.data.frame(x)
rownames(x) <- genes
colnames(x) <- c("H000", "H003", "H006",
                 "H009", "H012", "H018",
                 "H024", "H036", "H048",
                 "H072", "H120", "H168")
write.table(x, file = paste0(opt$outFolder, "/HMM.", opt$states, ".gene.matrix.tsv"),
            row.names = T, col.names = T, 
            quote=F, sep="\t")


# 5. save transition probabilities matrix
y <- matrix(getpars(fm1)[(nstates(fm1)+1):(nstates(fm1)^2+nstates(fm1))], byrow=TRUE,nrow=nstates(fm1),ncol=nstates(fm1))
write.table(y, file = paste0(opt$outFolder, "/HMM.", opt$states, ".transition.matrix.tsv"),
            row.names = T, col.names = T, 
            quote=F, sep="\t")


# 6. save response matrix
z <- as.data.frame(summary(fm1, which="response"))
write.table(z, file = paste0(opt$outFolder, "/HMM.", opt$states, ".response.matrix.tsv"),
            row.names = T, col.names = T, 
            quote=F, sep="\t")


# 7. save workspace
save.image(file = paste0(opt$outFolder, "/HMM.", opt$states, ".RData"), safe = TRUE)
