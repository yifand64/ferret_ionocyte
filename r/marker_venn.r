# given a few marker lists, generate a venn diagram showing genes that are in
# common. match based on ENSEMBL_ID and show gene symbol if possible.

cat(sprintf("Loading libraries.. \n"))
library(ggplot2)
library(VennDiagram)
# source("http://bioconductor.org/biocLite.R")
#biocLite()
# biocLite("limma")

args <- commandArgs(trailingOnly = TRUE)
n_args = length(args)
inputs = args[c(1:(n_args-1))]
output= args[n_args]
inputs_string = sprintf("%s", paste(inputs, collapse=" "))
cat(sprintf("Input datasets:       %s \n", inputs_string))
cat(sprintf("Output pdf:          %s \n", output))
pdf(output)

ensembl_cols = list()
gene_symbols = list()
dfs = list()

cat(sprintf("Loading data.. \n"))
for (i in 1:length(inputs)){
        m_list_file = inputs[i]
        cat(sprintf("Loading dataset #%i: %s.. \n", i,m_list_file))
        df <- read.delim(m_list_file)
        ensembl_cols[[i]] <- df[["ENSEMBL_ID"]]
        dfs[[i]] <- df
	head(ensembl_cols[[i]])

}

for (i in 1:length(ensembl_cols))
{
        cat(sprintf("%i genes in dataset %i \n",length(ensembl_cols[[i]]), i))
}

cat(sprintf("Counting intersections.. \n"))
#Draw the plot (assume only 2 datasets for the minute)

venn.plot <- draw.pairwise.venn(
area1 = length(ensembl_cols[[1]]),
area2 = length(ensembl_cols[[2]]),
cross.area = length(intersect(ensembl_cols[1], ensembl_cols[2])),
category = inputs,
fill = c("blue", "red"),
lty = "blank",
cex = 2,
cat.cex = 2,
cat.pos = c(285, 105),
cat.dist = 0.09,
cat.just = list(c(-1, -1), c(1, 1)),
ext.pos = 30,
ext.dist = -0.05,
ext.length = 0.85,
ext.line.lwd = 2,
ext.line.lty = "dashed") 
grid.draw(venn.plot)

