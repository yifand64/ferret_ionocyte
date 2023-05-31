#!/usr/bin/Rscript
#Run principal components analysis.
#Settings
take_log2		= TRUE
show_condition 	= FALSE
transpose_data 	= TRUE
title = "Single cell MiSeq \n marker genes"

homedir = path.expand("~")
source(paste(homedir, "/lib/CompoHeatMap/compoHeatMap.R", sep="/"))

cat(sprintf("Sourcing scripts..\n"))
scripts = c("samples.r",
			"pca.r",
			"qc.r")
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
this.script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dir <- dirname(this.script.name)

for (i in 1:length(scripts)) {
  	script.name = scripts[[i]]
  	other.name <- paste(sep="/", script.dir, script.name)
    print(paste("Sourcing", other.name, "from", this.script.name))
  	source(other.name)
    
}

cat(sprintf("Beginning PCA ..   \n"))

args <- commandArgs(trailingOnly = TRUE)
n_args = length(args)
if(n_args == 0)
{
	stop("ERROR: The first command line arg must be a counts table of expression data. ")
}
if(n_args > 1){
	#load qc data.
	qc_file = args[2]

}else{
	qc_file = NULL
}


# load counts for the genes of interest
expression_data_file = args[1]
#exclude = c("aneth", "ndo","ntero") # c("Paneth","Lgr5Hi","Endo")
exclude = c("ackground", "opulation")
expression_data = load_counts(expression_data_file, "normalised_count", exclude)#, drop_genes_file="cell_cycle_all.txt")
counts_data = expression_data[["counts"]] # load counts data for the samples 
complete_data = expression_data[["data"]]

if(!is.null(qc_file)){
	cat(sprintf("QC data is provided, loading it from: %s ..\n", qc_file))
	qc = load_qc(qc_file)
}else{
	qc = NULL
}

if(transpose_data){
	counts_data = t(counts_data)# load supplementary data for the samples (cells)
}

if(take_log2){
	print("Taking Log2+1  transform..")
	counts_data <- counts_data + 1
	counts_data <- log2(counts_data)
}


#cat(sprintf("Z-score normalising.. \n"))
#counts_data <- scale(counts_data)         # zscore -normalise
# clip = 3

# cat(sprintf("Clipping to Z = +- %s .. \n", clip))
# counts_data <- pmin(pmax(counts_data, -clip), clip)    
	
cat(sprintf("Defining cell type colours ..  \n"))
library(RColorBrewer)

# this should eventually be read from a file:
sample_type_labels = c("Paneth", "Lgr5Hi", "Lgr5lo", "Entero", "Endo")
sample_class_strings = sample_type_labels # strings must be detected in the sample name
sample_condition_strings = c("No.cell", "Population", "background","single")
sample_conditions = c("No cell background", "Population", "background","single cell")


## Define the colors corresponding to the cell types (or other types)
ctypes=sample_type_labels #paste0("type_", 1:4)
ct.colors=colorRampPalette(brewer.pal(9,"Set1"))(length(ctypes)+1) #c("#ffff33", "#ff9933", "#99ccff", "#33cc33")
names(ct.colors)=ctypes
## Sample colors
samples=sample_conditions
sample.colors=colorRampPalette(brewer.pal(9,"Set2"))(length(sample_conditions)+1)# c("grey", "black")
names(sample.colors)=samples

cat(sprintf("Labeling rows.. \n"))

## Define data for for row conditions => will be used to create a colorstack
# find sample classes:
cat(sprintf("Assigning samples to classes.. \n"))
sample_colors = colorRampPalette(brewer.pal(length(sample_type_labels), "Set1"))(length(sample_type_labels))
sample.groups = assign_samples(rownames(counts_data), 
									sample_type_labels, 
									sample_class_strings, 
									"Unknown", 
									sample_colors=sample_colors)

print("row conditions")
print(sample.groups)

#sample.groups=c(rep("type_1", nrow(counts_data))) #c(rep("type_1", 2), rep("type_2", 2), rep("type_3", 3), rep("type_4",1))
names(sample.groups) = rownames(counts_data)

if(show_condition){
	row.samples = assign_samples(rownames(counts_data), sample_conditions, sample_condition_strings, "Single cell")
	names(row.samples)=rownames(counts_data)
}

##PCA
info("Running PCA ..  \n")
pca = prcomp(counts_data)
#print(summary(pca))

info("Generating plots..")
visualise.pca(pca, sample.groups)

if(!is.null(qc))
{
	check.pcs(qc, pca, sample.groups, add.count.string=T)
}else{
	warn("No QC data, cannot check PCs for correlations with technical factors!")
}