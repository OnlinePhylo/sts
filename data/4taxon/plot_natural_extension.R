#!/usr/bin/env Rscript

library(ggplot2)
library(gridExtra)

# Theme
theme_set(theme_bw(20))
theme_update(legend.key=element_blank())

args <- commandArgs(TRUE)

stopifnot(length(args) == 3)

natural_ext <- read.csv(args[1], as.is=TRUE)
mb <- read.csv(args[2], as.is=TRUE, sep='\t', comment.char='[')
# Drop burn-in
mb <- tail(mb, nrow(natural_ext))
stopifnot(nrow(natural_ext) == nrow(mb))
comb <- cbind(natural_ext, mb)

to_burn <- floor(nrow(comb) * 0.25)
comb <- tail(comb, nrow(comb) - to_burn)

svg(args[3], height=6, width=15)
p1 <- ggplot(comb, aes(x=LnL, y=posterior, shape=leaves)) +
  geom_point(alpha=0.6) +
  xlab("MrBayes LnL") +
  ylab(expression(posterior["natural extension"]))
  #scale_color_gradient(name="Forest Length", low='red') +
p2 <- ggplot(comb, aes(x=TL, y=forest_length, color=LnL, shape=leaves)) +
  geom_point(alpha=0.6) +
  scale_color_gradient(name="MB LnL", low='red') +
  xlab("MrBayes total branch length") +
  ylab(expression(paste("Forest length on", R-1)))
grid.arrange(p1, p2, nrow=1)
dev.off()
