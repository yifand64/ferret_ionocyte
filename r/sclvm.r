# Run single cell latent variable models (scLVM)
# to estimate the effect of technical/biological/cell cycle
# variation, and to try to correct for variation
# caused by cell cycling.

# Written based on the tutorial / R markdown at:
# https://github.com/PMBio/scLVM/blob/master/R/tutorials/scLVM_vignette.Rmd



run_sclvm = function(norm.counts, 
						remove.cell.cycle=TRUE, 
						draw.plots=TRUE,
						cell_cycle_genes=NULL,
						write.corrected.output=TRUE,
						loess=FALSE)
{

	library(genefilter)
	library(statmod)
	library(scLVM)
	library(rPython)
	library(biomaRt)
	library(gplots)
	cat(sprintf("Starting scLVM.. \n"))
	cat(sprintf("Counts data:"))
	print(dim(norm.counts))
	# python.exec("from scLVM import scLVM")
	# python.exec("import scLVM as sclvm")
	cat(sprintf("Estimating technical noise.. \n"))

	pdf("scLVM_init.pdf", width=12, height=9)
	if(loess){
		tech.noise = fitTechnicalNoise(norm.counts, fit_type = 'logvar', use_ERCC = FALSE, plot=FALSE) 
	}else{
		tech.noise = fitTechnicalNoise(norm.counts, fit_type = 'log', use_ERCC = FALSE, plot=FALSE) 
	}
	cat(sprintf("Calling variable genes..\n"))
	is_het = getVariableGenes(norm.counts, tech.noise$fit, plot=TRUE)
	print(table(is_het))
  	cat(sprintf("Found %i variable genes.. \n", table(is_het)["TRUE"]))
	if(is.null(cell_cycle_genes))
	{
		#get cell cycle genes from GO 
		cat(sprintf("Getting cell cycle genes from GO.. \n"))
		ens_ids_cc <<- getEnsembl('GO:0007049')
    	# print(ens_ids_cc)
	}else{
	  ens_ids_cc = cell_cycle_genes
	}
  
  	#global variable definitions <<- fixes some 'object not found' errors.
	Y <<- t(log10(norm.counts+1)) #normalised trandformed read counts
	#Y <<- t(norm.counts)
	print("Set up scLVM counts matrix. Size: ")
	print(dim(Y))
	genes_het_bool <<- as.vector(is_het) #variable genes
	geneID <<- rownames(norm.counts) #gene IDs
  
	# mart <- useMart("ensembl")
	# listDatasets(mart=mart)
	# mart <- useDataset("mmusculus_gene_ensembl", mart = mart)
	# cat(sprintf("Converting to ENSEMBL ids..\n"))
	# genes.table <- getBM(filters= "ensembl_gene_id",
	#                      attributes= c("external_gene_id", "ensembl_gene_id", "description"), 
	#                       values= geneID, mart= mart)
	# print(genes.table)
  
	tech_noise = as.vector(tech.noise$techNoiseLog) #technical noise
	cat(sprintf("Building scLVM object..\n"))
	#construct and initialize new scLVM object
	sclvm = new("scLVM")
	cat(sprintf("Initialising..\n"))
	sclvm = scLVM::init(sclvm,Y=Y, tech_noise = tech_noise)
	cat(sprintf("Done.\n"))
	cat(sprintf("Fitting latent factors for %i cell-cycle genes..\n", length(ens_ids_cc)))
	CellCycleARD = scLVM::fitFactor(sclvm, 
	                          geneSet=ens_ids_cc, 
	                          k=20,
	                          use_ard = TRUE)
	cat(sprintf("Done.\n"))
	plot(seq(1, length(CellCycleARD$X_ard)), CellCycleARD$X_ard, xlab = '# Factor', ylab = 'Variance explained')
	title('Variance explained by latent factors')
	CellCycle = scLVM::fitFactor(sclvm, 
								geneSet = ens_ids_cc,
								k=1)
	#Get cell-cycle factor
	Kcc = CellCycle$K
	Xcc = CellCycle$X
	#Plot inferred similarity matrix
	image(Kcc,xaxt = "n", yaxt = "n", xlab = 'cells', ylab = 'cells')
	title('Similarity matrix based on cell cycle')
	idx_het = which(is_het)
	dev.off()


	# fit the model. uncomment idx=idx_het to fit for variable genes only.
	cat(sprintf("Running variance decomposition.. [may take a while]\n"))
	sclvm = varianceDecomposition(sclvm, K=Kcc) #, idx = idx_het)
  	cat(sprintf("Done. \n"))
	# get variance components
	results_var = getVarianceComponents(sclvm)
	var_filtered = results_var$var[results_var$conv,] # filter out genes for which vd has not converged
	head(var_filtered)

	# get corrected expression levels
	cat(sprintf("Getting corrected expression levels..\n"))
	Ycorr = getCorrectedExpression(sclvm)
	
	info("Saving corrected counts")
	# transform back to non-log counts
	exp10 = function(x){10^x}
	write.table(exp10(Ycorr) - 1, file="corrected_counts_scLVM.txt", sep="\t", quote=F)
	

	# if(write.corrected.output)
	# {
	# 	corr_counts_file = "corrected_counts.txt"
	# 	cat(sprintf("Writing corrected counts to %s ..\n", corr_counts_file))
	# 	corrected = data.frame(t(Ycorr))
	# 	print("Corrected counts df has size: ")
	# 	print(dim(corrected))
	# 	print(head(corrected))
	# 	rownames(corrected) = rownames(norm.counts)
	# 	colnames(corrected) = colnames(norm.counts)
	# 	print("After renaming columns and rows:")
	# 	print(dim(corrected))
	# 	print(head(corrected))
	# 	write.table(corrected, file=corr_counts_file, sep="\t", quote=F)
	# }
	#dim(Ycorr)
	var_mean = apply(var_filtered,2,mean)
	print(var_mean)

	if(!draw.plots){
		
	}else{
		
		pdf("scLVM_plots.pdf")
		colors = c('Green','Blue','Gray')
		pie(var_mean, , col = colors)
		idx_lmm = idx_het[1:5]
		
	  	cat(sprintf("Fitting LMM without correction..\n"))
		# fit lmm without correction
		res_nocorr = LMM(sclvm, K = NULL,idx = idx_lmm,verbose=TRUE)
		
		# fit lmm with correction
		cat(sprintf("Fitting LMM with correction..\n"))
	      
		res_corr = LMM(sclvm, K = Kcc, idx = idx_lmm,verbose = TRUE)
	  
	  	cat(sprintf("Generating heatmaps.. \n"))
		heatmap.2(res_nocorr$beta, Rowv = NULL, Colv = NULL, dendrogram = "none",
		          labCol = as.character(idx_lmm), labRow = as.character(idx_lmm),srtCol = 0, key=T,density.info = "none",
		          trace="none", breaks=seq.int(from = -0.6, to = 1.0, length.out = 13), main = 'Without Correction')
		
		heatmap.2(res_corr$beta, Rowv = NULL, Colv = NULL, dendrogram = "none",
		          labCol = as.character(idx_lmm), labRow = as.character(idx_lmm),srtCol = 0, key=T,density.info = "none",
		          trace="none", breaks=seq.int(from = -0.6, to = 1.0, length.out = 13), main = 'With Correction')
	  
		Yhet = Y[,idx_het]
		geneSymbols = getSymbols(colnames(Yhet))
		
		
		gene_plot = "Gata3"
		idx_gene = which(geneSymbols==gene_plot)
		cat(sprintf("Running PCA on corrected data.. \n"))
		#PCA on corrected data
		pcaCorr = prcomp(Ycorr,2)

		cat(sprintf("Drawing plots..\n"))
		print(Ycorr[,idx_gene])

		print(pcaCorr$x[,1])

		d <- qplot(pcaCorr$x[,1], pcaCorr$x[,2], xlab = 'PC1', ylab = 'PC2') + ggtitle('PCA corrected gene expression') 
		#d <- qplot(pcaCorr$x[,1], pcaCorr$x[,2], colour=as.factor(Ycorr[,idx_gene]), xlab = 'PC1', ylab = 'PC2') +ggtitle('PCA corrected gene expression') + scale_color_continuous(name =gene_plot)
		print(d)

		cat(sprintf("Running PCA on uncorrected data.. \n"))
		#PCA on uncorrected data
		pca = prcomp(Yhet,2)
		#d <- qplot(pca$x[,1], pca$x[,2], colour=Yhet[,idx_gene], xlab = 'PC1', ylab = 'PC2')+ggtitle('PCA uncorrected gene expression') + scale_color_continuous(name =gene_plot)
		d <- qplot(pca$x[,1], pca$x[,2], xlab = 'PC1', ylab = 'PC2')+ggtitle('PCA uncorrected gene expression') 
		print(d)

	}

	
	return (exp10(Ycorr) - 1)
}