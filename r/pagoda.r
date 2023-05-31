


do.pagoda <- function(raw.counts, knn.error.models, n.cores=1, use.msigdb=T, use.only=0, test=F)
{
	knn = knn.error.models
	cd = raw.counts
	#normalize out expected levels of technical and intrinsic biological noise
	info("Normalising variance")
	varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = n.cores, plot = TRUE)
	info("Top overdispersed genes")
	print(sort(varinfo$arv, decreasing = TRUE)[1:10])
	# Controlling for sequencing depth
	info("Controlling for number of detected genes")
	varinfo <- pagoda.subtract.aspect(varinfo, colSums(cd[, rownames(knn)]>0))

	info("Loading gene sets")
	if(!use.msigdb){
		library(org.Mm.eg.db)
		# translate gene names to ids
		info("Loading GO pathways")
		ids <- unlist(lapply(mget(rownames(cd), org.Mm.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
		info(sprintf("Mapped %s gene names to ids (%s genes discarded)", length(ids[!is.na(ids)])))
		rids <- names(ids); names(rids) <- ids
		# convert GO lists from ids to gene names
		
		if(use.only > 0){
			gos.interest <- unique(c(ls(org.Mm.egGO2ALLEGS)[1:use.only]))
		}else{
			gos.interest <- unique(ls(org.Mm.egGO2ALLEGS))
		}
		
		go.env <- lapply(mget(gos.interest, org.Mm.egGO2ALLEGS), function(x) as.character(na.omit(rids[x])))
		go.env <- clean.gos(go.env) # remove GOs with too few or too many genes
		go.env <- list2env(go.env) # convert to an environment
		set.type = "Gene Ontology"
	}else{
		
		if(test){
			warn("Test mode is ON. not using the complete MSigDB genesets, only like 100")
			load("MSigDB_test.Rdata")
			go.env = msig.test
			set.type = "MSigDB-Test-Only"
		}else{
			info("Loading mSigDB (from Livnat and Matan)")
			load("~/ref/gene_lists/MSigDB_mm10.RData")
			go.env = msig.Mm
			# puts MSigDB gene sets into a (R enivronment) called 'go.env'
			set.type = "MSigDB"
		}
	}
	n.sets = length(go.env)

	info(sprintf("Running pagoda weighted pca on %s gene sets [%s]", n.sets, set.type))
	pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components=1, n.cores=n.cores)
	return(pwpca)

	# info("Clustering genes to find variable 'de novo' gene sets")
	# clpca <- pagoda.gene.clusters(varinfo, trim = 7.1/ncol(varinfo$mat), n.clusters=50, n.cores=n.cores, plot = TRUE)

	# get full info on the top aspects
	tam <- pagoda.top.aspects(pwpca, n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))
	# determine overall cell clustering
	hc <- pagoda.cluster.cells(tam, varinfo)

	pagoda.show.pathways(c("GO:0022008","GO:0048699"), varinfo, go.env, 
		cell.clustering = hc, margins = c(1,5), 
		show.cell.dendrogram = TRUE, showRowLabels = TRUE, showPC = TRUE)



	# return(pwpca)
	
	# df <- pagoda.top.aspects(pwpca, return.table = TRUE, plot = TRUE, z.score = 1.96)
	# head(df)



	
	
	# info("Finding significantly variable gene sets")
	# df <- pagoda.top.aspects(pwpca, clpca, return.table = TRUE, plot = TRUE, z.score = 1.96)

	#clustering:
	# tam <- pagoda.top.aspects(pwpca, clpca, n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))
	# # determine overall cell clustering
	# hc <- pagoda.cluster.cells(tam, varinfo)


	return(tam)
}