# simple script, load a complete table with counts, FKPM, TPM etc, and spit out only
# the normalised_count

cat(sprintf("Starting ..   \n"))

args <- commandArgs(trailingOnly = TRUE)
n_args = length(args)
if(n_args == 0)
{
	stop("ERROR: The first command line arg must be a counts table of expression data. ")
}



expression_data_file = args[1]
exclude = c("ackground", "opulation")
count_string = "normalised_count"

cat(sprintf("Loading expression data from  --> %s\n", expression_data_file))
data <- read.delim(expression_data_file)
colnames = colnames(data)
samples = colnames[ grep(paste(count_string,collapse="|"), colnames)]
samples = samples[ grep(paste(exclude, collapse="|"), invert=TRUE, samples)]
cat(sprintf("Found %i samples: \n", length(samples)))
cat(sprintf("		%s \n", samples))
x = strsplit(samples, count_string)
keeps <- c("GENE_SYMBOL","transcript_id.s.","ENSEMBL_ID", samples)

data = data[keeps]
data = data[, c(1,2,3,7,8,9,6,5,4)]
write.table(data, file=paste(expression_data_file, "_trimmed.txt", sep=""), sep="\t",  row.names = F)