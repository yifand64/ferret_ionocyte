



cell.cycle.phase <- function(s.obj, draw.plots=T)
{
	# source("~/dev/adam/rna_seq/r/util.r")
	# script("seurat")
	# library(plyr)
	all = s.obj
	
	info("Running unsupervised clustering of cell cycle stages")
	x = fetch.data(all, c("G1.S_mean_normalised", "G2.M_mean_normalised"))
	all@data.info$Cell.Cycle.Stage.Unsup = pam(x, k=4)$clustering
	all@data.info$Cell.Cycle.Stage.Unsup = mapvalues(all@data.info$Cell.Cycle.Stage.Unsup, 1:4, c("G0", "G1", "S", "G2"))
	info("Inferred cell-cycle stages:")
	print(table(all@data.info$Cell.Cycle.Stage.Unsup))
	if(draw.plots)
	{
		info("Drawing plot [new cols]")
		scatter(all, "G1.S_mean_normalised", "G2.M_mean_normalised", colour.by="Cell.Cycle.Stage.Unsup", 
			pdf="cc_unsup.pdf", width=11, height=9, cols=perceptual_rainbow_16[c(1, 5, 10, 15)])
	}
	
	ncells = length(all@cell.names)
	info("Running supervised clustering")
	# thresholds are a list of x/y thresholds that define regions in the G1/G2 space
	thresh=list("G0"=c(0,0), "G1"=c(0,0), "S"=c(0,0))
	#thresh=list("G0"=c(0,0), "G1"=c(0,-0.3), "S"=c(0.5,-0.3))

	# print(thresh$G0[1])
	# print(sum(all@data.info$G1.S_mean_normalised<-1))
	g0 = which(all@data.info$G1.S_mean_normalised<thresh$G0[1] & all@data.info$G2.M_mean_normalised<thresh$G0[2])
	#print(g0)
	#return(FALSE)
	g1 = which(all@data.info$G1.S_mean_normalised>thresh$G1[1] & all@data.info$G2.M_mean_normalised<thresh$G1[2])

	s = which(all@data.info$G1.S_mean_normalised>thresh$S[1] & all@data.info$G2.M_mean_normalised>thresh$S[2])

	g2.m = which(!(1:ncells) %in% c(g0, g1, s))

	cc = rep("Unknown", ncells)
	cc[g0] = "G0"
	cc[g1] = "G1"
	cc[s] = "S"
	cc[g2.m] = "G2.M"
	info("Inferred cell-cycle stages:")
	all@data.info$Cell.Cycle.Stage.Sup = cc
	if(draw.plots){
		info("Drawing plot [new cols]")
		p = scatter(all, "G1.S_mean_normalised", "G2.M_mean_normalised", colour.by="Cell.Cycle.Stage.Sup",
			 cols=perceptual_rainbow_16[c(1, 5, 10, 15)]) + geom_vline(xintercept=0, linetype="dotted") + geom_hline(yintercept=0, linetype="dotted") 
		ggsave(p, filename="cc_sup.pdf", width=11, height=9,)
	}
	return(all)
}

# draw a bar plot that shows the fraction of each cluster classified
# into each phase of the cell cycle
cell.cycle.state.by.group <- function(s.obj, group.by, cex.all=16, xlab="", ylab="Proportion of cells")
{
	ggplot(s.obj@data.info, aes_string(x=group.by, fill="Cell.Cycle.Stage.Sup")) + geom_bar(position = "fill") + 
	theme_bw() + scale_fill_manual("CC-stage", values=perceptual_rainbow_16[c(1, 5, 10, 15)]) + 
	theme(axis.line.x = element_line(color="gray15", size = 0.75, lineend="round"),
   			axis.line.y = element_line(color="gray15", size = 0.75, lineend="round")) + 
    theme(
        panel.border= element_blank(),
        axis.title.x = element_text(vjust=-0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size=cex.all, colour="gray22")) + xlab(xlab) + ylab(ylab)
}

#assumes that cell.cycle.phase was called first
order.cells <- function(s.obj, draw.plots=T)
{
	info("Ordering cells")
	
	x = fetch.data(s.obj, c("Cell.Cycle.Stage.Sup", "G1.S_mean_normalised", "G2.M_mean_normalised"))
	print(table(x$Cell.Cycle.Stage.Sup))
	#order cells in G0 by expression of G1
	g0 = x[x$Cell.Cycle.Stage.Sup == "G0",]
	g0 = g0[order(g0$G1.S_mean_normalised),]
	g0$rank = 1:nrow(g0)
	
	print(sprintf("G0 cells ordered from %s", paste(range(g0$rank), collapse="-")))
	#print(head(g0, n=20))

	#order cells in G1 by expression of G1
	g1 = x[x$Cell.Cycle.Stage.Sup == "G1",]
	g1 = g1[order(g1$G1.S_mean_normalised),]
	g1$rank = 1:nrow(g1) + nrow(g0)
	print(sprintf("G1 cells ordered from %s", paste(range(g1$rank), collapse="-")))
	
	#order cells in S by expression of G2
	s = x[x$Cell.Cycle.Stage.Sup == "S",]
	s = s[order(s$G2.M_mean_normalised),]
	s$rank = 1:nrow(s) + (nrow(g0) + nrow(g1))
	print(sprintf("S cells ordered from %s", paste(range(s$rank), collapse="-")))

	#print(s)

	#order cells in G2 inversely by expression of G2
	g2 = x[x$Cell.Cycle.Stage.Sup == "G2.M",]
	g2 = g2[order(g2$G2.M_mean_normalised, decreasing=T),]
	#print(g2)
	g2$rank = 1:nrow(g2) + (nrow(g0) + nrow(g1) + nrow(s))
	print(sprintf("G2 cells ordered from %s", paste(range(g2$rank), collapse="-")))


	rval = rbind.data.frame(g0, g1, s, g2)
	rval = rval[match(s.obj@cell.names,rownames(rval)),]

	s.obj@data.info$Cell.Cycle.Order = rval$rank
	if(draw.plots){
		info("Drawing plot")
		scatter(s.obj, "G1.S_mean_normalised", "G2.M_mean_normalised", colour.by="Cell.Cycle.Order", cols=perceptual_rainbow_16,
			pdf="cc_ordered.pdf", width=11, height=9)

	}
	return(s.obj)
}


# get a list of genes, ranked by their association (pearson's correlation)
# with CC
cc.genes <- function(s.obj, cc.geneset)
{
	if(is.null(cc.geneset))
	{
		error("Cannot run test for CC-associated genes without a set of CC genes! Provide a cc.geneset")
		return(FALSE)
	}
	cc.score = unlist(fetch.data(s.obj, "Cell.Cycle_mean_normalised"))
	complexity.score = unlist(fetch.data(s.obj, "nGene"))
	#print(cc.score)
	cc.cor = apply(s.obj@data, 1, function(x){cor(unlist(x), cc.score)})
	complexity.cor = apply(s.obj@data, 1, function(x){cor(unlist(x), complexity.score)})
	d = data.frame(cc.cor, complexity.cor)
	d$gene = rownames(d)
	d$Mean.Expr = rowMeans(s.obj@data)
	#d$assign = pam(d, k=3)$clustering
	return(d)
}



# Evans (1996) suggests for the absolute value of r:
# ===================================================
# .00-.19 “very weak”  
# .20-.39 “weak” 
# .40-.59 “moderate” 
# .60-.79 “strong” 
# .80-1.0 “very strong”

# Identify cell-cycle associated PCs
check.pcs <- function(
					s.obj,
					cc.geneset=NULL,
					use.n.genes.hypergeometric=100,
					p.threshold=0.001,
					corr.threshold=0.5,
					pca=NULL, 
					groups=NULL,
					n_components=20)
{
	
	initial_wd = getwd()
	if(is.null(cc.geneset))
	{
		error("Cannot run test for enriched genes without a set of CC genes! Set use.score=F or provide a cc.geneset")
		return(FALSE)
	}
	check_dir = "Cell_Cycle_Associations"
	
	dir.create(check_dir, showWarnings = FALSE)
	setwd(check_dir)

	if(is.null(pca)){
		info("Running PCA"); 
		pca = prcomp(t(s.obj@data), scale=F, center=T)
	}
	if(is.null(groups)){groups=factor(rep("All cells", nrow(pca$x)))}

	scores <- data.frame(groups, pca$x[,1:n_components])
	colnames(scores)[1] = "groups"
	scores$sample_name <- rownames(scores) 

	# print(head(scores))

	cat(sprintf("Checking correlation of PCs with Cell-cycle: \n"))
	# cat(sprintf("QC data: \n"))
	# print(head(qc))
	metrics = c("Cell.Cycle_mean_normalised", "G1.S_mean_normalised", "G2.M_mean_normalised", "MHC.II_mean_normalised")
	metric_names = c("Cell.Cycle_mean_normalised"="Cell-cycle score", "G1.S_mean_normalised"="G1/S phase score", "G2.M_mean_normalised"="G2/M phase score", "MHC.II_mean_normalised"="MHCII")
	
	qc = fetch.data(s.obj, metrics)
	qc$sample_name = s.obj@cell.names

	qc = merge(qc, scores, by="sample_name")
	
	ok = vector()
	lg = get.loaded.genes(pca, n_genes=use.n.genes.hypergeometric, components=1:n_components)
	#print(head(lg))
	p = vector()

	m = NULL

	for(j in 1:n_components)
	{
		ok[j] = T
		pearson_cor = vector()
		for(i in 1:length(metrics))
		{
			qc_metric = metrics[[i]]
			pc = paste("PC",j, sep="")
			
			#cat(sprintf("Checking correlation of %s with %s .. \n", qc_metric, pc))
			pearson = round(cor(qc[qc_metric], qc[pc], method="pearson"), 3)
			spearman = round(cor(qc[qc_metric], qc[pc], method="spearman"), 3)
			#cat(sprintf("Spearman = %f, Pearson = %f \n", spearman, pearson))
			ncols = length(unique(groups))
			pearson_cor[i] = pearson
			if(abs(spearman) > corr.threshold | abs(pearson) > corr.threshold)
			{
				ok[j] = F
				warn(sprintf("%s is associated with %s! Plotting.", pc, qc_metric))
				correlation_label = paste("Pearson = ", pearson,"\n", "Spearman = ", spearman ,sep="")
				pl = ggplot(qc, aes_string(x=qc_metric, y=pc, colour=groups)) + geom_point() + theme_bw()  +
					scale_colour_manual("Cell type",  values=default.cols(ncols))  + ggtitle(correlation_label) + xlab(metric_names[qc_metric])
				ggsave(pl, file=paste("Correlation", qc_metric,pc, ".pdf",sep="_"), height=6, width=6)
			}	
		}

		if(is.null(m))
		{
			cors = pearson_cor
			names(cors) = unlist(metrics)
			cors = data.frame(t(cors))

		}else{
			cors = rbind(cors, pearson_cor)
		}
		#print()

		# Hypergeometric test of enrichment of Cell-cycle genes in the loadings for a given PC.
		# https://www.biostars.org/p/15548/
		# phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
		# where q=size of overlap-1; m=number of upregulated genes in experiment #1; 
		# n=(total number of genes on platform-m); k=number of upregulated genes in experiment #2.
		
		n = nrow(pca$rotation)
		cc.geneset = cc.geneset[!is.na(cc.geneset)]
		m = length(cc.geneset)

		loaded.genes = unlist(unname(as.character(lg[, j])))
		#print(loaded.genes)
		olap = toupper(as.character(cc.geneset)) %in% toupper(loaded.genes)
		olap.genes = cc.geneset[olap]
		olap.str = paste(olap.genes, collapse=", ")
		no = length(olap.genes)	
		if(no>0)
		{
			info(sprintf("%s CC genes were found loaded highly on PC-%s [%s]", no, j, olap.str))
		}
		
		q = sum(olap) -1
		k = length(loaded.genes)
		p[j]  = phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
		info(sprintf("Hypergeometric test [PC-%s]. N=%s, overlap-1=%s, m=%s, k=%s. Pval=%s", j, n, q,m, k, signif(p[j], 3)))
	}

	
	info("Adjusting p-values (FDR)")
	p.fdr = p.adjust(p, method="fdr")
	p.bonf = p.adjust(p, method="bonferroni")
	ok.genes = p.bonf > p.threshold
	ok.score = apply(cors, 1, function(x){any(abs(x)>corr.threshold)})

	rval = data.frame(1:n_components, p.fdr, p.bonf, !ok.genes, ok.score)
	colnames(rval) = c("PC", "p.fdr", "p.bonf", "CC-genes-loaded", "CC-Correlates")
	rval = cbind.data.frame(rval, cors)
	
	rval = cbind.data.frame(rval, t(summary(pca)$importance[, 1:n_components]))

	info("PC --  Cell Cycle  check complete")

	cc.associated = which(!ok)
	usable = which(ok)
	info(sprintf("PCs %s are cell-cycle associated and PCs %s are not", paste(cc.associated, collapse=", "), paste(usable, collapse=", ")))
	
	setwd(initial_wd)
	write.table(rval, file="cell_cycle_assocations.txt", sep="\t", quote=F, row.names=F)
	return(rval)
}