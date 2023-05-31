
# library(Matrix)
# library(Rcpp)
# library(bit64)
library(data.table)
library(rhdf5)
# library(Rtsne,Rtsne)
# library(dplyr,mutate)
# library(irlba,irlba)


# must be run in a folder containing three files.
# 1. counts_table.txt -- a cells/barcode matrix in tab-delimited text
# 2. barcodes.txt  -- a tab-delimited list of barcodes, the column names for the matrix
# 3. genes.txt -- a tab-delimited list of gene names, the row names
name_10x_table <- function()
{
	library(useful)
	cat("Loading unnamed table \n")
	cts = read.delim("counts_table.txt", check.names=F, header=F)
	
	cat("Loading gene names (from genes.txt) \n")
	genes = read.delim("genes.txt", header=F)

	cat("Loading cell names from barcodes.txt \n")
	b = read.delim("barcodes.txt", check.names=F, header=F)

	cat("Naming table \n")
	colnames(cts) = as.character(unlist(b))
	rownames(cts) = make.names(as.character(unlist(genes)), unique = T)

	cat("Finished naming table. Writing\n")
	#print(corner(cts))
	write.table(cts, file="counts_table_named.txt", sep="\t", quote=F)
}

merge.list <- function(x, by="GENE_SYMBOL")
{
	info("Doing merge.list")
	rval = x[[1]]
	n = length(x)
	for(i in 2:n)
	{
		info(sprintf("Merge %s of %s, current size: %s", i-1, n-1, paste(dim(rval), collapse="x")))
		info(sprintf("i=%s", i))
		if(is.null(dim(x[[i]]))){
			warn(sprintf("Table #%s is empty! Skipping it", i))
		}else{
			print(dim(x[[i]]))
			new.cols = grep(by, colnames(x[[i]]), invert=T)
			if(nrow(rval) == nrow(x[[i]]))
			{
				#info("Here")
				
				#print(new.cols)
				info(sprintf("	adding %s new columns", length(new.cols)))
				rval = cbind(rval, x[[i]][, new.cols], by=by)
			}else{
				stop("Merge (rbind) failed! Tables have differing number of rows!")
				

				#rval = merge(rval, x[[i]], by=by)
			}
		}
	}
	info("Merge (rbind) complete!")
	return(rval)
}


# load a set of 10X experiments from different counts tables. readr::read_tsv is faster
# than read.delim, and write_tsv / write_rds is faster than write.table.
merge.batch <- function(sample.file, use.readr=T, progress=T, to.rds=T)
{
	library(useful)
	library(parallel)
	if(use.readr){library(readr)}

	sf = read.delim(sample.file, check.names=F, stringsAsFactors=F)
	files = as.character(sf[,2])
	batch.names = unlist(as.character(sf[,1]))
	source("~/dev/adam/rna_seq/r/samples.r")
	info("Checking provided file paths")
	if(!(check.paths(files)))
	{
		stop("Exiting. Fix samples file")
	}

	info(sprintf("Loading %s 10X batches:", length(batch.names)))
	library(parallel)
	cl <- makeCluster(8L)
	clusterExport(cl, c("info", "flog.info")) #, "iono", "gam"))
	print(batch.names)
	if(progress){
		library(pbapply)
		counts = pblapply(files, read.ct, cl = cl, fast=use.readr)
	}else{
		counts = mclapply(files, read.ct, mc.cores = detectCores(), fast=use.readr)
	}
	

	info(sprintf("Loaded %s counts tables", length(counts)))
	# grab gene_symbols and ensembl_ids from the first counts table.
	info("Got gene ids:")
	gene_ids = counts[[1]][,1:2]

	# there are some genes like, IL11RA2 that occur in multiple places in the
	# genome, and therefore have the same symbol, and different ENSEMBL ids. so for
	# these, we paste on the ENSEMBL id to resolve ambiguity (also R won't allow duplicate
	# rownames)
	
	same_gene = which(duplicated(gene_ids$GENE_SYMBOL))
	info(sprintf("Resolving %s duplicate ids by pasting on ENSEMBL ids", length(same_gene)))
	gene_ids$use = gene_ids$GENE_SYMBOL
	resolved.ids = paste(as.character(gene_ids$use)[same_gene], 
									as.character(gene_ids$ENSEMBL)[same_gene], sep="_")
	print(resolved.ids)

	gene_ids$use[same_gene] = resolved.ids
	# print("Added:")
	# print(grep("_", as.character(gene_ids$use), value=T))

	use.ids = as.character(gene_ids$use)
	still.same = use.ids[duplicated(use.ids)]

	info(sprintf(" %s duplicate ids remain", length(still.same)))
	if(length(still.same)>0){
		print(still.same)
	}

	info(sprintf("Loaded %s count tables", length(counts)))

	info("Building groups vector")
	n.samples.in.each.dataset = as.numeric(unlist(lapply(counts, ncol)))	
	groups = factor(unlist(mapply(rep, batch.names, n.samples.in.each.dataset))) 
	print(table(groups))

	info("Merging tables..")
	merged = merge.list(counts)
	
	# this make.names call was a bad idea, added X's to the front of many (mostly RIK) gene names
	#rownames(merged) = make.names(as.character(unlist(merged$GENE_SYMBOL)), unique = T) 
	
	rownames(merged) = gene_ids$use
	merged$Row.names = NULL


	colnames(merged) = paste(colnames(merged), groups, sep="_")

	# remove some unnecessary columns the merge added

	# print(grep("by", colnames(merged), value=T))
	# print(grep("GENE", colnames(merged), value=T))
	# print(grep("ENSEMBL", colnames(merged), value=T))
	#write.table(merged, file="merged_counts_dirty.txt", sep="\t", quote=F)

	merged[,grep("by", colnames(merged), value=T)] <- NULL
	merged[,grep("GENE", colnames(merged))] <- NULL
	merged[,grep("ENSEMBL", colnames(merged))] <- NULL

	info(sprintf("Built %s merged table. Table contains %s missing values (NAs)", paste(dim(merged), collapse="x"), sum(is.na(merged))))
	info("Saving..")
	if(use.readr){
		if(to.rds){
			write_rds(merged, "merged_counts.rds", compress="gz")
		}else{
			readr::write_tsv(merged, path="merged_counts.txt", col_names=T)
			write.table(rownames(merged), file="merged_counts_gene_names.txt", quote=F, row.names=F)
		}
	}else{
		if(to.rds)
		{
			saveRDS(merged, file="merged_counts.rds")
		}else{
			write.table(merged, file="merged_counts.txt", sep="\t", quote=F)
		}
		
	}
	info("Done!")
	return(list(merged=merged, groups=groups))
}

# read a counts table
read.ct <- function(f, fast=T)
{
	info(sprintf("Reading from %s", f))
	if(fast){ct = readr::read_tsv(f, progress=F)}else{
		ct = read.delim(f, check.names=F, stringsAsFactors=F)
	}
	
	# info("converting to dataframe")
	# ct = data.frame(ct)
	# if("X" %in% colnames(ct))
	# {
	# 	warn("Replacing screwed up gene column name")
	# 	genes = unlist(as.character(ct$X))
	# 	rownames(ct) = make.names(genes, unique = T)
	# 	ct$X = NULL
	# }
	#colnames(ct)[1] <- "GENE_SYMBOL"
	#print(corner(ct))
	ct

}





# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

decompress_seqs_cpp <- function(binary_seq, seq_len) {
    .Call('cellrangerRkit_decompress_seqs_cpp', PACKAGE = 'cellrangerRkit', binary_seq, seq_len)
}


#' Load molecule info h5 file from the Cell Ranger pipeline
#'
#' @param gbm GeneBCMatrix object to load into
#' @param barcode_len Length of cell barcode in nucleotides
#' @import rhdf5
#' @import bit64
#' @import Rcpp
#' @import data.table
#' @export
#' @return A GeneBCMatrix object containing the loaded molecule info
get_molecule_info_complete <- function(gbm, barcode_len=14) {
  h5_path <- file.path(gbm@pipestance_path, "outs", "molecule_info.h5")
  required_cols <- c("barcode", "gem_group", "gene", "umi", "reads", "nonconf_mapped_reads", "unmapped_reads") #"umi_corrected_reads", "barcode_corrected_reads", "conf_mapped_uniq_read_pos",  
  	
  cols <- h5ls(h5_path)$name
  for (col in required_cols) {
    if (!(col %in% cols)) {
      stop(sprintf("Missing '%s' column in molecule info h5", col))
    }
  }
  col_data <- lapply(required_cols, function(col) {
    as.integer(h5read(h5_path, col, bit64conversion='bit64'))
  })
  names(col_data) <- required_cols
  col_data[["barcode"]] <- decompress_seqs_cpp(as.integer(col_data[["barcode"]]),
                                               seq_len=barcode_len)  # only works for up to 15mers

  
  cd = data.frame(col_data, stringsAsFactors=F)
  # return(cd)
  # print(colnames(cd))
  # print(head(cd))
  # print(range(cd$unmapped_reads))
  # print(colSums(cd[1:100000, c("reads", "unmapped_reads")]))
  dt <- as.data.table(cd)
  
  # Aggregate reads per-cell (barcode) vcfr5rtr5
  setkey(dt, barcode, gem_group, gene, umi)
  dt <- dt[, j=list(reads=sum(reads), unmapped_reads=sum(unmapped_reads), nonconf_mapped_reads=sum(nonconf_mapped_reads)), by=list(barcode, gem_group, gene, umi)]
  #return (dt)

  # # Aggregate reads per-molecule
  # 
  # dt <- dt[, j=list(reads=sum(reads)), by=c('barcode', 'gem_group', 'gene', 'umi')]
  gbm@molecule_info <- dt
  return(gbm)
}



### note that the $reads field is the number of transcriptome-mapped reads.
### see cellranger description: http://software.10xgenomics.com/single-cell/pipelines/latest/output/molecule_info
### "Number of reads that confidently mapped to this putative molecule."
mapping.rates <- function(gbm)
{
	dt = gbm@molecule_info
	dt <- dt[, j=list(reads=sum(reads), unmapped_reads=sum(unmapped_reads), nonconf_mapped_reads=sum(nonconf_mapped_reads),
		ngene=length(unique(gene)), numi=length(unique(umi))), by=list(barcode)]

	total = sum(dt$reads) / (sum(dt$reads) + sum(dt$unmapped_reads) + sum(dt$nonconf_mapped_reads))

	info(sprintf("Total: %s ", total))
	pass.filter.bcs = gsub("-1", "", colnames(gbm@mat))
	dt$trans.mapping.rate = dt$reads / (dt$unmapped_reads + dt$reads + dt$nonconf_mapped_reads)
	dt$genome.mapping.rate = (dt$reads + dt$nonconf_mapped_reads) / (dt$unmapped_reads + dt$reads + dt$nonconf_mapped_reads)
	dt = dt[dt$reads>0 | dt$unmapped_reads>0,]
	info(sprintf("Transcriptome mapping. Mean: %s and median: %s (all cells)", mean(dt$trans.mapping.rate), median(dt$trans.mapping.rate)))
	info(sprintf("Genome mapping. Mean: %s and median: %s (all cells)", mean(dt$genome.mapping.rate), median(dt$genome.mapping.rate)))
	dt.bcpass = dt[dt$barcode %in% pass.filter.bcs, ]


	info(sprintf("Mean Reads per Cell: %s ", mean(dt.bcpass$reads + dt.bcpass$unmapped_reads + dt.bcpass$nonconf_mapped_reads)))
	info(sprintf("Median Genes per Cell: %s ", median(dt.bcpass$ngene)))
	info(sprintf("Median UMI per Cell: %s ", median(dt.bcpass$numi)))

	dt$pass.filter = dt$barcode %in% pass.filter.bcs

	info(sprintf("Mean: %s and median: %s (pass-filter cells)", mean(dt.bcpass$trans.mapping.rate), median(dt.bcpass$trans.mapping.rate)))
	
	d = as.data.frame(dt)

	d = d[d$reads > 100,] # dont plot really shitty cells
	min.genes=800
	min.trans.mapped=0
	g = ggplot(d, aes(x=ngene, y=trans.mapping.rate)) + geom_vline(xintercept=min.genes) + theme_bw() + 
					geom_hline(yintercept = min.trans.mapped) + xlab("Number of genes detected") + ylab("% of reads mapping to transcriptome") + 
					geom_polygon(aes(fill = ..level..), stat = "density2d", alpha = .2) + 
					scale_fill_continuous("Density of Cells", low = "#56B1F7", high = "#132B43", guide=FALSE ) + 
					geom_point(size=0.5, aes(colour=pass.filter)) + 
					scale_colour_manual("Pass CellRanger filter", values=default.cols(2)) + 
					theme(axis.line.x = element_line(color="gray15", size = 0.75, lineend="round"),
						axis.line.y = element_line(color="gray15", size = 0.75, lineend="round"),
					    panel.border = element_blank(),
					    panel.background = element_blank(), 
					    strip.background = element_blank(), 
					    strip.text.x = element_text(size = 16, face="bold"),
					    text = element_text(size=16, colour="gray22"))
	print(g)

	return(dt)
}


### gbm -- gene barcode matrix, a cellranger class with everything in it.
visualise.gene <- function(gene_name, gbm)
{
	gi = gbm@molecule_info[gbm@molecule_info$gene==grep(gene_name, gbm@gene_symbols),]

	n.umis = length(unique(gi$umi))
	n.reads = sum(gi$reads)
	n.cells = length(unique(gi$barcode))
	info(sprintf("Gene %s is detected in %s cells. There are %s UMIs and %s reads.", gene_name, n.cells, n.umis, n.reads))

	dt <- gi[, j=list(reads=sum(reads), umis=length(unique(umi))),by=list(barcode)]
	df = as.data.frame(dt)
	df$reads.log10 = log10(df$reads+1)
	df$umis.log10 = log10(df$umis+1)

	#print(head(df))
	dfm = melt(df[, c("reads.log10", "umis.log10", "barcode")])
	g = ggplot(dfm, aes(x=value, fill=variable)) + geom_density(alpha=0.5) + xlim(c(0,4)) + theme_bw() + 
		ggtitle(sprintf("%s expression \n Detected in %s cells out of %s cells (10X)", gene_name, n.cells, length(unique(gbm@molecule_info$barcode))))
	print(range(df$reads.log10))
	print(g)
	return(df)
}

rank.genes.by.reads <- function(gbm)
{
	dt = gbm@molecule_info
	dt <- dt[, j=list(reads=sum(reads)), by=list(gene)]
	df = as.data.frame(dt)
	df = df[order(df$reads, decreasing=T),]
	df$Gene_symbol = gbm@gene_symbols[df$gene]
	return(df)
}

### gbm -- gene barcode matrix, a cellranger class with everything in it.
setup_tenx <- function(path_to_pipestance, genome="mm10")
{
	library(cellrangerRkit)
	#analysis <- load_cellranger_analysis_results(path_to_pipestance)
	gbm <- load_cellranger_matrix(path_to_pipestance, barcode_filtered=T, genome = genome)
	gbm <- get_molecule_info_complete(gbm)
	return(gbm)
	#genes = read.delim(paste0(path_to_pipestance, "/outs/filtered_gene_bc_matrices/mm10/genes.tsv"), header=F)
	#colnames(genes) = c("ENSEMBL", "Gene_symbol")
	#return(list(gbm=gbm, genes=genes))
}





