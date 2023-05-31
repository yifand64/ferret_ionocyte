# run Monocle to search for differentiation state tree structure in single-cell
# rna-seq data.

library(monocle)
library(reshape)

 #settings:
count_string				= "normalised_count"
expr_in_n_cells 			= 50


# args <- commandArgs(trailingOnly = TRUE)
# n_args = length(args)
# if(n_args == 0)
# {
# 	stop("ERROR: The first command line arg must be a counts table of expression data. ")
# }

# initial.options <- commandArgs(trailingOnly = FALSE)
# file.arg.name <- "--file="
# script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
# script.basename <- dirname(script.name)

# # source other scripts in the same directory:
# scripts = c("samples.r")
# for (i in 1:length(scripts)) {
# 	script_name = scripts[[i]]
# 	other.name <- paste(sep="/", script.basename, script_name)
# 	source(other.name)
# 	print(paste("Sourcing",other.name,"from",script.name))
# }


# # load counts for the genes of interest
# expression_data_file = args[1]
# exclude = c("ackground", "opulation","Paneth", "Endo") # c("Paneth","Lgr5Hi","Endo")

# expression_data = load_counts(expression_data_file, "normalised_count", exclude)
# counts_matrix = expression_data[["counts"]] # load counts data for the samples 
# complete_data = expression_data[["data"]]

# cat(sprintf("Loaded data matrix. Dimensions: \n"))
# print (dim(counts_matrix))



setup.monocle <- function(seurat.obj)
{
	pd <- new("AnnotatedDataFrame", data = seurat.obj@data.info)
	counts_matrix = seurat.obj@data
	print(pd)

	#no annotation of genes for the moment:
	#fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
	cd <- newCellDataSet(as.matrix(counts_matrix), phenoData = pd, expressionFamily=negbinomial.size())

}



monocle <- function(seurat.obj, markers=c("Lgr5", "Hspd1","Ldha","Mttp","Adck3","Ces2c","Sis","Cdc20",
											"Agr2","Arl6ip1","Clca3","Ccnb2","Tpsg1","Defa-rs1","Defa-rs7","Gm10104",
											"Fbxo22","Hck","Tnfrsf19","Kctd12","Clca3","Pygl","Agr2","Ccl6"))
{

	cat(sprintf("Starting Monocle pipeline.. \n"))

	#an AnnotatedDataFrame object, where rows are cells, and columns are cell attributes (such as cell
	#type, culture condition, day captured, etc.)
	# cat(sprintf("Building sample sheet frame: \n"))
	# sample_class_title = "Cell type"
	# sample_class_labels = c("Paneth","Lgr5Hi","Lgr5lo","Entero", "Endo") #c("Organoid", "In vivo") #"CD103+ CD11b+", "CD103+ CD11b-", "CD103- CD11b+", "CD103- CD11b-", "Macrophage"
	# sample_class_strings = sample_class_labels	# c("Org", "SI")  #"_CD103posCD11bpos_DC", "_CD103pos_","_CD11bpos_","_DN", "macro")
	# #sample_colors = colorRampPalette(brewer.pal(length(sample_class_labels), sample_color_palette))(length(sample_class_labels))
	# sample_assignments = assign_samples(colnames(counts_matrix), sample_class_labels, sample_class_strings, "Unknown")
	# sample_sheet <- data.frame(sample_assignments)
	# colnames(sample_sheet) = c("FACs_type")
	# rownames(sample_sheet) = colnames(counts_matrix)
	# print(sample_sheet)

	pd <- new("AnnotatedDataFrame", data = seurat.obj@data.info)
	counts_matrix = seurat.obj@data
	print(pd)

	#no annotation of genes for the moment:
	#fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
	cd <- newCellDataSet(as.matrix(counts_matrix), phenoData = pd)

	cat(sprintf("Constructed CellDataSet: \n"))
	print(cd)


	cat(sprintf("Detecting the expression of genes expressed in at least %i cells in the dataset.. \n", expr_in_n_cells))
	cd <- detectGenes(cd, min_expr = 0.1)
	expressed_genes <- row.names(subset(fData(cd), num_cells_expressed >= expr_in_n_cells))
	#print(expressed_genes)

	cat(sprintf("Checking quality control: \n"))
	#valid_cells <- row.names(subset(pData(cd), Cells.in.Well == 1 & Control == FALSE & Clump == FALSE & Debris == FALSE & Mapped.Fragments > 1000000))
	#print(table(valid_cells))
	#cd <- cd[,valid_cells]
	#print(head(pData(cd)))

	# Log-transform each value in the expression matrix.
	#L <- log(exprs(cd[expressed_genes,]))

	L <- exprs(cd[expressed_genes,])

	# Standardize each gene, so that they are all on the same scale,
	# Then melt the data with plyr so we can plot it easily"
	melted_dens_df <- melt(t(scale(t(L))))
	# Plot the distribution of the standardized gene expression values.
	pdf("distribution_of_gene_expression.pdf",height=8.5,width=11)
	qplot(value, geom="density", data=melted_dens_df) + stat_function(fun = dnorm, size=0.5, color='red') +
	xlab("Standardized log(FPKM)") +
	ylab("Density") + theme_bw()
	dev.off()

	# eventually here we can search for differentially expressed genes, or load in a list of markers. for now just look at 
	# genes that we think are important: [top 3 SCDE per cluster]
	marker_genes <- row.names(subset(fData(cd), rownames(cd) %in% markers))

	# cat(sprintf("Running differential expression test..\n"))
	# diff_test_res <- differentialGeneTest(cd[marker_genes,])


	#plot expression of a couple of interesting genes:
	interesting = c("Lgr5", "Clca3", "Defa-rs1")
	pdf("jitter.pdf", width=11, height=8.5)
	jitter <- cd[row.names(subset(fData(cd),
						rownames(cd) %in% interesting)),]
	plot_genes_jitter(jitter, grouping="orig.ident", ncol=length(interesting))
	dev.off()

	# for the moment keep this simple, only use marker genes that are expressed in expr_in_n_cells 
	# to try to build the ordering tree.
	cat(sprintf("Using ordering genes: \n"))
	ordering_genes <- intersect(marker_genes, expressed_genes)
	cd <- setOrderingFilter(cd, ordering_genes)
	print(ordering_genes)

	cat(sprintf("Running ICA.. \n"))
	cd <- reduceDimension(cd, use_irlba=FALSE)

	cat(sprintf("Running Monocle.. \n"))
	cd <- orderCells(cd, num_paths=2, reverse=TRUE)
	
	cat(sprintf("Plotting spanning tree.. \n"))
	pdf("spanning_tree.pdf", width=11, 8.5)
	plot_spanning_tree(cd, show_cell_names=TRUE)
	dev.off()


	cat(sprintf("Plotting gene expression dynamics in pseudotime.. \n"))
	cd_filtered <- cd[expressed_genes, pData(cd)$State != 3]
	my_genes <- row.names(subset(fData(cd_filtered), rownames(cd_filtered) %in% c("Lgr5", "Hspd1", "Ldha")))
	cds_subset <- cd_filtered[my_genes,]
	pdf("genes_in_psuedotime.pdf", width=11, height=8.5)
	plot_genes_in_pseudotime(cds_subset, color_by="orig.ident")
	dev.off()

}













