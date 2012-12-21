#!/usr/bin/env Rscript

library(ggplot2)
library(gridExtra)

# Theme
theme_set(theme_bw(20))
theme_update(legend.key=element_blank())


args <- commandArgs(TRUE)

stopifnot(length(args) == 3)

cherries <- read.csv(args[1], as.is=TRUE)
mb <- read.csv(args[2], as.is=TRUE, sep='\t', comment.char='[')
stopifnot(nrow(cherries) == nrow(mb))
comb <- cbind(cherries, mb)

to_burn <- floor(nrow(comb) * 0.25)
comb <- tail(comb, nrow(comb) - to_burn)

svg(args[3], height=6, width=15)
p1 <- ggplot(comb, aes(x=LnL, y=posterior, color=forest_length)) +
  geom_point(alpha=0.6) +
  scale_colour_gradient(name="Forest Length", low='red') +
  xlab("MrBayes LnL") +
  ylab(expression(posterior[cherries]))
p2 <- ggplot(comb, aes(x=TL, y=forest_length, color=LnL)) +
  geom_point(alpha=0.6) +
  scale_colour_gradient(name="MB LnL", low='red') +
  xlab("MrBayes total branch length") +
  ylab("Forest length on cherries")
grid.arrange(p1, p2, nrow=1)
dev.off()
