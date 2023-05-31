#
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
	cat(sprintf("Loading %s.. \n", m_list_file))
        df <- read.delim(m_list_file)
	ensembl_cols[[i]] <- df[["ENSEMBL_ID"]]
	dfs[[i]] <- df
	
}

for (i in 1:length(ensembl_cols))
{
	cat(sprintf("%i genes in dataset %i \n",length(ensembl_cols[[i]]), i))
}

cat(sprintf("Counting intersections.. \n"))
universe <- Reduce(union, ensembl_cols)
print(universe)
Counts <- matrix(0, nrow=length(universe), ncol=length(inputs))
colnames(Counts) <- inputs
for (i in 1:length(universe))
{
 	cat(sprintf("Searching for %s .. \n", universe[i]))
	for(j in 1:length(inputs))
	{
		Counts[i,j] <- universe[i] %in% ensembl_cols[j]
		
	}
}
print(Counts)

vennDiagram(vennCounts(Counts)
