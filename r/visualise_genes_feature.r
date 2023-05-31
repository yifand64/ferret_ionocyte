# given a list of files, for each file,
# load a list of genes, and display their expression
# patterns in a single cell data set using seurat's feature plots.

draw_feature_plots(markers.file, 
					seurat.obj, 
					output_file=NULL)
{
	cat(sprintf("Reading markers from %s ..\n", markers.file))
	markers = read.delim(markers.file)
	genes = markers$GENE_SYMBOL
	pdf(output_file)
	feature.plot()
}