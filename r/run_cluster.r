# Run hierarchical clustering 

# Written by Adam, Jan 25, 2015.
#
# R script to cluster (probably rna-seq) expression data 
# Info on clustering and heatmap.2: https://www.biostars.org/p/18211/


library(methods)
	
# should switch to named args:
# https://cwcode.wordpress.com/2013/04/16/the-joys-of-rscript/

args <- commandArgs(trailingOnly = TRUE)
n_args = length(args)
if(n_args == 0)
{
	stop("ERROR: The first command line arg must be a counts table of expression data. ")
}
if(n_args > 1){
	title = args[2]
}else{
	title = "HiSeq data " # "Gut DCs and Macrophages"
}
if(n_args > 2){
	output_file=args[3]
}
if(n_args > 3){
	qc_file = args[4]
}else{
	qc_file = NULL
}

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

# source other scripts in the same directory:
scripts = c("samples.r",
	"cluster.r")
for (i in 1:length(scripts)) {
	script_name = scripts[[i]]
	other.name <- paste(sep="/", script.basename, script_name)
	source(other.name)
	print(paste("Sourcing",other.name,"from",script.name))
}

#Settings
#-------------------------------------------------------------------
# include NMF:
#script.dir <- dirname(sys.frame(1)$ofile)
#source(paste(script_dir,"nmf.r",sep="/"))

# settings:
NMF 					= FALSE 	# not sure if this is complete, definitely not tested.
n_de 					= 100
find_de_genes 			= FALSE 	# find the n_de most differentially expressed genes using edgeR
correlation_matrix 		= TRUE
tmm_normalise 			= FALSE
do_pca 					= FALSE


#---- which samples to include?
count_string = "normalised_count" #"FPKM")		#"expected_count")#,"counts")
exclude = c("ackground","opulation")
expression_data_file = args[1]

if(length(args) > 3){
	complete_expression_data_file = args[4]
	cat(sprintf("Loading the complete counts data frame from %s \n", complete_expression_data_file))
	complete_counts_data=load_counts(complete_expression_data_file, count_string, exclude)[["counts"]]
}

expression_data = load_counts(expression_data_file, count_string, exclude)
cat(sprintf("Counts loaded. \n"))
counts = expression_data[["counts"]] # load counts data for the samples 
complete_data = expression_data[["data"]]

# comment in to manually reorder columns:
#counts_data = counts_data[, rev(c(1,5,2,4,3))]

if("Type"  %in% colnames(complete_data))
{
	cat(sprintf("Can use \'Type\' information to annotate genes..\n"))
}else{
	cat(sprintf("WARN: No \'Type\' column. Cannot label gene types.. \n"))
	show_gene_type = FALSE
}
sample_names = colnames(counts)
# check for nans:
num_nans = sum(is.nan(as.matrix(counts)))
cat(sprintf("Counts data frame contains %i NaN values.. \n", num_nans))

if(do_pca)
{
	cat(sprintf("Running PCA ..\n"))
	pdf("pca_plots.pdf")
	pca <- prcomp(counts, scale=T)
	print(pca)
	cat(sprintf("Generating PCA plots.. \n"))

}

#correlation matrix
if(correlation_matrix){
	pdf("correlation_matrix.pdf")
	cat(sprintf("Generating correlation matrix ..\n"))
	cor_mat = melt(cor(counts, method="pearson"))
	rownames(cor_mat) = NULL
	cor_mat_plot = ggplot(cor_mat, aes(X1, X2, fill = value)) + geom_tile() + scale_fill_gradientn(colours=rev(colorRampPalette(brewer.pal(9,"YlGnBu"))(50    ))) 
	cor_mat_plot = cor_mat_plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
	print(cor_mat_plot)
}

# differentially expressed genes:
cat(sprintf("\"Finding differentially expressed genes is %s ...\"\n", if(find_de_genes) "ON" else "OFF"))
if(find_de_genes){
	#refactored this and haven't tested it:
	find_differential_genes_edger(counts, sample_names, n_genes=n_de)
	
}

if(tmm_normalise){
	counts = normalise.edgr.tmm(counts)
}
# can show column sums to see differences in library size:
#colSums(y$pseudo.counts)

#run non-negative matrix factorization to discover covarying modules.
if(NMF){
	print("Running non-negative matrix factorization..")
	run_nmf(counts, 3, "nmf.pdf")
}
# print(head(counts))
cat(sprintf("Calling cluster function..\n"))
cluster(counts, 
	show_gene_type = show_gene_type,
	complete_data=complete_data, 
	pdf.output=TRUE)

