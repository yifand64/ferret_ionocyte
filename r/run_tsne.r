

# settings:
count_string				= "normalised_count"
pca_reduce_dimensionality 	= FALSE
scde.correct 				= TRUE
num_pcs 					= 8


args <- commandArgs(trailingOnly = TRUE)
n_args = length(args)
if(n_args == 0)
{
	stop("ERROR: The first command line arg must be a counts table of expression data. ")
}

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

# source other scripts in the same directory:
scripts = c("samples.r", "tsne.r","scde.r")
for (i in 1:length(scripts)) {
	script_name = scripts[[i]]
	other.name <- paste(sep="/", script.basename, script_name)
	source(other.name)
	print(paste("Sourcing",other.name,"from",script.name))
}

# load counts for the genes of interest
expression_data_file = args[1]
exclude = c("ackground", "opulation","Paneth_A60","Paneth_A72") # c("Paneth","Lgr5Hi","Endo")

expression_data = load_counts(expression_data_file, "expected_count", exclude)
counts = expression_data[["counts"]] # load counts data for the samples 
complete_data = expression_data[["data"]]

cat(sprintf("Loaded data matrix. Dimensions: \n"))
print (dim(counts))

cat(sprintf("Starting tSNE pipeline.. [Convert to PCA scores is %s] \n", if(pca_reduce_dimensionality) "ON" else "OFF"))



if(pca_reduce_dimensionality){
	counts = log2(counts+1)
	tsne_input = convert_to_pca_scores(counts, pcs=num_pcs)
	cat(sprintf("After PCA dimensionality reduction, data matrix: \n"))
	print (dim(tsne_input))
}else{
	if(scde.correct)
	{
		# keep.only = 20
		# cat(sprintf("Only keeping %i cells ..\n",keep.only))
		# counts = counts[, sample(ncol(counts), keep.only)]

		dist.type = "reciprocal"
		scde.models = scde.compute.error.models(counts)
		scde.plot.dropout(scde.models[["error.models"]], scde.models[["expression.prior"]])

		if(identical(dist.type, "reciprocal"))
		{
			tsne_input = t(scde.dist.reciprocal(scde.models[["counts"]], scde.models[["error.models"]], scde.models[["expression.prior"]]))
		}else
		{
			if(identical(dist.type, "mode.rel"))
			{
				tsne_input = t(scde.dist.mode.rel(scde.models[["counts"]], scde.models[["error.models"]], scde.models[["expression.prior"]]))
			}else
			{
				if(identical(dist.type, "direct"))
				{
					tsne_input = t(scde.dist.direct.dropout(scde.models[["counts"]], scde.models[["error.models"]], scde.models[["expression.prior"]]))
				}
			}
		}

		print("Size of scde-corrected distance matrix:")
		print(dim(tsne_input))
		print(head(tsne_input))
	}else
	{
		tsne_input = counts
	}
	

}

do.tsne(tsne_input, n=tsne_run_n_times, parallel=T)

 





