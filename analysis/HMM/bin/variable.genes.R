# 1. set working directory
setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/")


# 2. read expression matrix of 12448 genes
expression.matrix <- read.table("expression/QN.merged/selected.genes.rep.2.3.after.QN.merged.tsv",
                                h=T, sep="\t")

# 3. retrieve set of not expressed genes
not.expressed.genes <- read.table("expression/silent.genes.txt", h=F, sep="\t", 
                                  stringsAsFactors = F)
not.expressed.genes <- not.expressed.genes$V1


# 4. retrieve set of DE genes
DE.genes <- read.table("expression/QN.merged/expression.matrix.tsv", h=T, sep="\t")
DE.genes <- rownames(DE.genes)


# 5. retrieve set of stably expressed genes 
stably.expressed.genes <- setdiff(rownames(expression.matrix), 
                                  c(not.expressed.genes, DE.genes))

# 6. create a df with gene label (not expressed, stably expressed, DE)
x <- data.frame(gene_id = c(DE.genes, stably.expressed.genes, not.expressed.genes),
                group = c(rep("DE", length(DE.genes)), 
                          rep("stably expressed", length(stably.expressed.genes)),
                          rep("not expressed", length(not.expressed.genes))))

# 6. change to HMM wd

df <- data.frame(stringsAsFactors = F)

runs <- c("marks.expression", "marks")

for ( k in runs ) {
  
  setwd(paste0("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/HMM/",
        k))
  
  for ( i in 2:20 ) {
    
    # 7.11. read matrix of states' sequence for each gene
    m <- read.table(paste0("HMM.", i, ".gene.matrix.tsv"), 
                    h=T, sep="\t")
    
    # 7.12. count, for each gene, how many states are found
    m$type <- apply(m, 1, function(x){length(table(x))})
    m$type <- ifelse(m$type == 1, "stable", "variable")
    
    m$gene_id <- rownames(m)
    m <- merge(m, x, by = "gene_id")
    
    y <- melt(table(m[, c("type", "group")]))
    y$n_states <- i
    y$total <- c(rep(8030, 2), rep(1552, 2), rep(2666, 2))
    y$run <- k
    
    df <- rbind(df, y)
    
  }
  
}




df$group <- factor(df$group, levels = c("not expressed", 
                                        "stably expressed",
                                        "DE"))
df$run <- gsub("marks.expression", "marks&expression", df$run)


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/HMM/variable.genes.pdf",
    width = 12, height = 4)
ggplot(df[df$type == "variable", ],
       aes(x=n_states, y = value/total, group = run,
           color = run)) +
  geom_point() +
  geom_line() +
  facet_grid(~group) +
  scale_y_continuous(labels = percent_format(),
                     limits = c(0, 1)) +
  ylab("% of variable genes") +
  xlab("# of states")  +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        strip.text.x = element_text(size = 14),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        strip.background.x = element_blank()) +
  scale_x_continuous(breaks = 2:20)

dev.off()
  