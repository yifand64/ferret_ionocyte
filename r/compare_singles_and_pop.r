3# Written by Adam, Jan 25, 2015.
#
# R script to cluster rna-seq expression data for dataset provided by
# Omer via Yarden. Dataset is read counts (normalised by Yarden), then
# filtered according to:
# 
# Info on clustering and heatmap.2: https://www.biostars.org/p/18211/


# include NMF:
#script.dir <- dirname(sys.frame(1)$ofile)
#source(paste(script_dir,"nmf.r",sep="/"))

#Settings
#-------------------------------------------------------------------
NMF=FALSE

args <- commandArgs(trailingOnly = TRUE)
n_args = length(args)

if(n_args == 0)
{
	stop("ERROR: The first command line arg should be a counts table of expression data. ")

}

# heatmap title:
if(n_args > 1){
	title = args[2]
}else{
	title = "HiSeq data " # "Gut DCs and Macrophages"
}

#load other scripts
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

# source other scripts in the same directory:
scripts = c("samples.r","examine_subpopulation.r")
for (i in 1:length(scripts)) {
	script_name = scripts[[i]]
	other.name <- paste(sep="/", script.basename, script_name)
	source(other.name)
	print(paste("Sourcing",other.name,"from",script.name))
}

library(reshape)
library(edgeR)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(data.table)

# settings
tmm_normalise = TRUE
take_log2 = TRUE

cell_types = c("Entero", "Endo","Paneth") #Lgr5Hi", "Lgr5Lo") #

cat(sprintf("Comparing single cells and population averages for  \n"))
cat(sprintf("%s  \n", cell_types))


#---- which samples to include?
count_string = c("normalised_count") #"FPKM")		#"expected_count")#,"counts")
exclude = c("ackground")

expression_data = load_counts(args[1], count_string, exclude)
cat(sprintf("Counts loaded. \n"))
counts_data = data.frame(expression_data[["counts"]]) # load counts data for the samples 
complete_data = expression_data[["data"]]


if(length(args) > 1)
{
	population_data_provided = TRUE
	population_data_file = args[2]
	cat(sprintf("Loading population data from %s \n", population_data_file))
	population_expression_data = load_counts(population_data_file, count_string, exclude)
	cat(sprintf("Population data loaded. \n"))
	population_counts_data = data.frame(population_expression_data[["counts"]]) # load counts data for the samples 
	population_complete_data = population_expression_data[["data"]]
	
}else
{
	population_data_provided = FALSE
	
}

generate.plot <- function(compare, single.av.col, pop.col, single_colnames=NULL)
{
	if(!is.null(single_colnames))
	{
		cat(sprintf("Single cell data columns: \n"))
		cat(sprintf("-----------------------------\n"))
		cat(sprintf(" %s  \n", single_colnames))
	}
	
	if(!population_data_provided){
		cat(sprintf("population_col: %s .. \n", population_colname))
	}
	cat(sprintf("single_cell_average_col: %s .. \n", single_cell_average_colname))
	cat(sprintf("-----------------------------\n \n \n"))

	pearson_cor = cor(compare[single_cell_average_colname], compare[population_colname], method="pearson")
	spearman_cor = cor(compare[single_cell_average_colname], compare[population_colname], method="spearman")

	cat(sprintf(" 		Correlations: Pearson= %f, Spearman= %f .. \n", pearson_cor, spearman_cor))
	cat(sprintf("		Generating scatter plot .. [%s vs %s]", single_cell_average_colname, population_colname))
	cat(sprintf("\n \n"))
	xlab = paste("Single cell average - ", cell_type, " Log2(counts + 1)")
	if(population_data_provided){
		 ylab = paste("Population [HiSeq] - ", cell_type, " Log2(counts + 1)")
	}else{
		 ylab = "Population Control Log2(counts + 1)"
	}
	
	xmin = 0
	xmax = 15
	ymin = 0
	ymax = 15

	correlation_label = paste("Pearson = ", pearson_cor,"\n", "Spearman = ", spearman_cor ,sep="")
	print("Single cell average:")
	print(head(compare[single_cell_average_colname]))
	print("population")
	print(head(compare[population_colname]))

	p = ggplot(compare, aes_string(x=single_cell_average_colname, y=population_colname)) + geom_point(size=0.75) + theme_bw()  +
	scale_x_continuous(xlab, breaks =  seq(xmin, xmax, 2), limits=c(xmin, xmax)) + 
	scale_y_continuous(ylab, breaks =  seq(ymin, ymax, 2), limits=c(ymin, ymax)) 
	p = p + annotate(geom="text", x=12, y=1, label=correlation_label, size=3, family="Times")
	print(p)
}

for (i in 1:length(cell_types)) {
		cell_type = cell_types[i]
		single_cell_cols = grep(paste("normalised_count_", cell_type,"_B", sep=""), colnames(counts_data))
		# print(single_cell_cols)
		cat(sprintf("Comparing %s [%i cells].. \n", cell_type, length(single_cell_cols)))
		
		pdf(paste(cell_type,".pdf", sep=""))
		
		# calculate the single cell average
		single_cell_average_colname = paste("normalised_count_", cell_type,"_singles_average", sep="")
		single_colnames = colnames(counts_data)[single_cell_cols]
		compare = complete_data[, c("GENE_SYMBOL","ENSEMBL_ID")]
		cat(sprintf("Adding single cell average column.. \n"))
		compare[single_cell_average_colname] = rowMeans(counts_data[single_cell_cols])
		print(head(compare))
		cat(sprintf("Done. \n"))

		# put the population data into the 'compare' dataframe
		cat(sprintf("Looking for population data.. \n"))
		if(population_data_provided){
			cat(sprintf("Population data provided.. (second arg)\n"))
			
			# # hack because of different sample naming in hiseq-miseq data sets!
			# if(cell_type=="Lgr5Hi")
			# {
			# 	cell_type="Lgr5_Hi"
			# }else if(cell_type=="Lgr5lo")
			# {
			# 	cell_type="Lgr5_lo"
			# }
			cat(sprintf("CellType= %s \n", cell_type))


			population_col_index = grep(paste("normalised_count_SI_", cell_type, sep=""), colnames(population_counts_data))
			population_colname = colnames(population_counts_data)[population_col_index]

			cat(sprintf("Using population data in column: %s \n", population_colname))
			cat(sprintf("Adding IDs.. \n"))

			population_counts_data["GENE_SYMBOL"] = population_complete_data["GENE_SYMBOL"]
			population_counts_data["ENSEMBL_ID"] = population_complete_data["ENSEMBL_ID"]

			cat(sprintf("Converting to data.table.. \n"))
			data_table_1 = data.table(compare)
			data_table_2 = data.table(population_counts_data)
			cat(sprintf("Merging.. \n"))
			# print("Table 1:")
			# print(head(data_table_1))
			# print("Table 2:")
			# print(head(data_table_2))
			cat(sprintf("Converting to back to data.frame.. \n"))
			# find the set of genes that are expressed in both datasets.
			compare = data.frame(merge(data_table_1, data_table_2, by="GENE_SYMBOL")) 
		}else
		{
			cat(sprintf("No population data provided, will search for population controls in the single cell dataset (first arg)..\n"))
			population_col_index = grep(paste("normalised_count_", cell_type,"_Pop", sep=""), colnames(counts_data))
			if(length(population_col_index > 1))
			{
				pop_i = 1
				cat(sprintf("Multiple population columns found!\n"))
				print(population_col_index)
				cat(sprintf( "Using %i ..\n", pop_i))
				
				population_col_index = population_col_index[pop_i]
			}
			population_colname = colnames(counts_data)[population_col_index]
			compare[population_colname] = counts_data[population_colname]
		}

		sample_columns = colnames(compare)[ grep(paste(count_string,collapse="|"), colnames(compare))]
		if(take_log2){
			print("Taking Log2+1  transform..")
			compare[sample_columns] <- compare[sample_columns] + 1
			compare[sample_columns] <- log2(compare[sample_columns])
		}

		if(tmm_normalise){
			print("TMM normalizing..")
			y2 = DGEList(counts=compare[sample_columns])
			y2 = calcNormFactors(y2)
			print(y2$samples) # show normalization factors:
			y2 = estimateCommonDisp(y2, verbose=TRUE)
			compare[sample_columns] = y2$pseudo.counts
		}

		# print(head(counts_data[,single_cell_average_colname], n=10))
		# print(head(counts_data[,population_colname], n=10))

		generate.plot(compare, single.av.col = single_cell_average_colname, pop.col=population_colname)
		

		cat(sprintf("|-----------------| 	|------------------------|		|---------------------|\n \n \n \n"))

}



