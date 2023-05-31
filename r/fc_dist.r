#install.packages("directlabels", repos = "http://r-forge.r-project.org")
library(directlabels)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
cat(sprintf("Input dataset:       %s \n", args[1]))
cat(sprintf("Output pdf: 	  %s \n", args[2]))
pdf(args[2])

asNumeric <- function(x) as.numeric(levels(x))[x]#as.character(x))
factorsNumeric <- function(d) modifyList(d, lapply(d[, sapply(d, is.factor)],
                                                   asNumeric))

fold_change_dist <- read.delim(args[1])

str(fold_change_dist)
fold_change_dist[,c('FC.Threshold')] = as.numeric(as.character(fold_change_dist[,c('FC.Threshold')]))
str(fold_change_dist)

number_ticks <- function(n) {function(limits) pretty(limits, n)}
p <- ggplot(fold_change_dist, aes(x=FC.Threshold, y=Number.of.genes, label=Samples, colour=Samples)) + geom_line() + scale_y_log10(breaks=number_ticks(10)) +scale_x_log10() + theme_bw() + theme(legend.text=element_text(size=0.4)) + geom_dl(aes(label=Samples, colour=Samples), list('last.points', cex = 0.4, hjust = 1))
print(p)

#direct.label(p, list(last.points, hjust = 0.7, vjust = 1))


