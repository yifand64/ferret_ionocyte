
# needed for normalization:
#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
library(edgeR)
library(gplots)
library(RColorBrewer)
library(fastcluster)
#library(vegan)
#library(amap)
#library(cairo):
library(ggplot2)
library(scales)

cat(sprintf("Sourcingx scripts..\n"))
scripts = c("samples.r")
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
this.script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dir <- dirname(this.script.name)
initial_wd = getwd()


for (i in 1:length(scripts)) {
  	script.name = scripts[[i]]
  	other.name <- paste(sep="/", script.dir, script.name)
    print(paste("Sourcing", other.name, "from", this.script.name))
  	source(other.name)
}


# the first and only arg is the expression data set.
args <- commandArgs(trailingOnly = TRUE)
n_args = length(args)
count_string = c("normalised_count","counts")
exclude = c("bulk", "Bulk")
expression_data_file = args[1]

# find sample names:
#data <- read.delim(expression_data_file)

#if there are columns to be ignored: 
# samples = colnames[ grep(paste(count_string,collapse="|"), colnames)]
# samples = samples[ grep(paste(exclude, collapse="|"), invert=TRUE, samples)]
# #cat(sprintf("\"Found column names  --> %s\"\n", colnames))
# keeps <- c("GENE_SYMBOL", samples)
# cat(sprintf("\"Selecting column:  %s\"\n", keeps))
# norm_data <- data[keeps]



count_string = "normalised_count"
exclude = c("ackground","opulation")
              
# Load RSEM counts table
expression_data = load_counts(expression_data_file, 
                              count_string, 
                              exclude, 
                              trim.sample.names=F)
              
info("Counts loaded.")
norm_data = expression_data[["counts"]]
# TMM normalization:
#print("TMM normalizing..")
#y <- DGEList(counts= norm_data[data_cols])
#y <- calcNormFactors(y)
#y = estimateCommonDisp(y, verbose=TRUE)
#print(y$samples)
#norm_data[data_cols] = y$pseudo.counts
#colSums(y$pseudo.counts)

# samples = colnames(norm_data)
# data_cols <- samples
# # log2(counts+1)

# info("Taking Log transform..")
# log_data = norm_data
# log_data[data_cols] <- log2( norm_data[data_cols] + 1)

# #print(str( log_data))
# info("Building data matrix")
# m <- as.matrix(log_data[data_cols])

# print(dim(norm_data))
# print(dim(log_data))
# print(head(log_data))

samples = colnames(norm_data)
log_data = data.frame(log2(norm_data + 1))
print(dim(norm_data))
print(dim(log_data))

fancy_scientific <- function(l) {
     # turn in to character string in scientific notation
     l <- format(l, scientific = TRUE)
     # quote the part before the exponent to keep all the digits
     l <- gsub("^(.*)e", "'\\1'e", l)
     # turn the 'e+' into plotmath format
     l <- gsub("e", "%*%10^", l)
     # return this as an expression
     parse(text=l)
}

info(sprintf("Generating histograms for %i samples..", length(samples)))
pdf(file="histograms.pdf")  
for (i in 1:length(samples))
{

		sample_name = samples[i]
		print(sample_name %in% colnames(log_data))

		
		info(sprintf("Sample  --> %s", sample_name))
		plot = ggplot(log_data, aes_string(x=sample_name)) + geom_histogram(aes(y=..density..)) + theme_bw() + xlab(paste("log_2 ",sample_name,sep="")) + scale_y_continuous(limits=c(0,2)) + geom_density(alpha=.2, fill="#FF6666")
		print(plot)
		#ggsave("histograms.pdf")
}

pdf(file="fold_change_dists.pdf") 
of_interest = paste("Entero",sep="|")
samples = samples[grep(of_interest, samples)]
pairs = as.matrix(expand.grid(samples, samples))
for (i in 1: nrow(pairs))
{
		pair = pairs[i,]
		sample_a = pair[1]
		sample_b = pair[2]
		if (sample_a != sample_b) {
			cat(sprintf("\"Generating fold change distribution for samples  %s vs %s\"\n", sample_a, sample_b))
			data_a = norm_data[sample_a]
			data_b = norm_data[sample_b]
			fc_vector = log2((data_a / data_b))
			d = data.frame(fc_vector)
			#print(d)
			plot = ggplot(d, aes_string(x=sample_a)) + geom_histogram(aes(y=..density..)) + theme_bw() + xlab(paste(sample_a,sample_b,sep=" vs ")) + geom_density(alpha=.2, fill="#FF6666") #+ scale_x_continuous(breaks = seq(-10,10,1))
			print(plot)
		}
		#ggsave("histograms.pdf")
}



