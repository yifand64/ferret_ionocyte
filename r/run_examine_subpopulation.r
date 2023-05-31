main <- function()
{
	library(ggplot2)
	library(RColorBrewer)
	# settings:
	count_string				  	= "normalised_count"
	population_of_interest_string  	= "No_cell"
	population_of_interest_cols	  	= list()
	size_of_random_sample		  	= 100
	n_diff_genes				 	= 100
	show_top_expressed_genes_also	= TRUE

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
	scripts = c("samples.r", "examine_subpopulation.r")
	for (i in 1:length(scripts)) {
		script_name = scripts[[i]]
		other.name <- paste(sep="/", script.basename, script_name)
		source(other.name)
		cat(sprintf(paste("Sourcing",other.name,"from",script.name, "\n")))
	}


	# load counts for the genes of interest
	expression_data_file = args[1]
	expression_data = load_counts(expression_data_file, "normalised_count")
	counts_matrix = expression_data[["counts"]] # load counts data for the samples 
	complete_data = expression_data[["data"]]

	cat(sprintf("Loaded data matrix. Dimensions: \n"))
	print (dim(counts_matrix))

	examine_subpopulation(counts_matrix, complete_data,
		#poi_cols=population_of_interest_cols, 
		name_of_population="No cell",
		poi_string = population_of_interest_string,
		n_background_samples = size_of_random_sample,
		n_diff_genes=n_diff_genes, 
		show_top=show_top_expressed_genes_also)


}

main()
