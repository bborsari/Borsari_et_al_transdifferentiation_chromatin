setwd("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/HMM/")

runs <- c("marks.expression", "marks")

x <- data.frame(stringsAsFactors = F)

for (k in runs) {
  
  tmp <- read.table(paste0("log.lik.HMM.", k))
  colnames(tmp) <- "logLik"
  tmp$run <- k
  tmp$n_states <- (1:nrow(tmp)) +1
  x <- rbind(x, tmp)
  
}

x$run <- gsub("marks.expression", "marks&expression", x$run)


pdf("~/public_html/Borsari_et_al_transdifferentiation_chromatin/HMM/logLik.pdf",
    width = 6, height = 4)
ggplot(x, aes(x=n_states, y=logLik, group=run, color=run)) +
  geom_point() +
  geom_line() +
  xlab("# of states")  +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        strip.text.x = element_text(size = 14),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        strip.background.x = element_blank()) +
  scale_x_continuous(breaks = 2:20)
dev.off()
