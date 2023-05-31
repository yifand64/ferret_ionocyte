library(Seurat)
library(reshape)
#sean simmons' function

# adam [march1,2017]: sean only returns the 'hurdle' component
# which is a combination (average) of the coefficient and p-values
# for the discrete and continuous parts of the model. fyi the discrete
# models the binary probability of a gene being 'on' (logistic regression)
# and continuous models its expression level as a guassian.


# note, x is output from mast. see de_test for specifics.
makeNice<-function(x,val="ko_vs_wtwt",comp="H")
{
	info("Cleaning up MAST output")
	#print(x)
	#x.use<-x[x$component==comp & x$contrast==val,]
	x.use <-x[x$component==comp & x$contrast==val,]
	info("here")
	pval<-x.use[,4]
	padj<-p.adjust(pval,"fdr")
	padj_strict<-p.adjust(pval,"BY")

	pval<-cbind(pval,padj)
	pval<-cbind(pval,padj_strict)


	rownames(pval)<-x.use[,1]
	x_logFC<-x[x$component=="logFC",]

	z<-cast(x_logFC,formula=primerid~component+contrast,value="z",fun.aggregate=sum)
	logfc<-cast(x_logFC,formula=primerid~component+contrast,value="coef",fun.aggregate=sum)
	rownames(z)<-z[,1]
	rownames(logfc)<-logfc[,1]
	cls<-colnames(z)
	for(i in 1:length(cls)){cls[i]=paste(cls[i],"_z",sep="")}
	colnames(z)<-cls

	z<-z[2:dim(z)[2]]
	logfc<-logfc[2:dim(logfc)[2]]

	#logfc<-cbind(logfc,z)

	lst<-x.use[,1]
	rownames(x.use)<-x.use[,1]
	#ord<-c()
	#num<-dim(z)[2]

	#for(i in 1:num){ord<-c(ord,i,(num+i))}

	#logfc<-logfc[,ord]

	ret<-cbind(pval,logfc)

	#ret<-ret[order(ret$pval),]

	ret["Gene"]<-rownames(ret)
	ret<-ret[c(dim(ret)[2],1:(dim(ret)[2]-1))]

	#print(head(ret))
	#adam added:
	colnames(ret) = gsub("pval", "p", colnames(ret))
	colnames(ret)[grep("logFC_g", colnames(ret))] <- "log2fc"
	###

	ret
}



mast_gsea <- function(zlm, group.var="groups", compare.group="pSC")
{
	library(GSEABase)


	#### get GSEA sets. this needs to be redone to be general.
	module <- "BTM"
	min_gene_in_module <- 5
	packageExt <- system.file("extdata", package='MAST')
	module_file <- list.files(packageExt, pattern = module, full.names = TRUE)
	gene_set <- getGmt(module_file)
	gene_ids <- geneIds(gene_set)
	gene_ids <- gene_ids[!names(gene_ids)%like%"TBA"&!names(gene_ids)%like%"B cell"]
	sets_indices <- limma::ids2indices(gene_ids, mcols(sca)$symbolid)
	# Only keep modules with at least min_gene_in_module
	sets_indices <- sets_indices[sapply(sets_indices, length) >= min_gene_in_module]


	gsea <- gseaAfterBoot(de@mast.glm, boots, sets_indices, hypothesis = CoefficientHypothesis(paste0(group.var, compare.group))) 
z_stat_comb <- summary(gsea, testType='normal')

	sigModules <- z_stat_comb[combined_adj<.01]
	gseaTable <- melt(sigModules[,.(set, disc_Z, cont_Z, combined_Z)], id.vars='set')
	ggplot(gseaTable, aes(y=set, x=variable, fill=value))+geom_raster() + scale_fill_distiller(palette="PiYG")	
}


