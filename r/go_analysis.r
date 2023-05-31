
source("gsea.r")
library(hash)
library(org.Mm.eg.db)
library(Mus.musculus)

jgc <- function() ##java garbage collection
{
	info("Garbage collecting..")
	gc()
	.jcall("java/lang/System", method = "gc")
} 

#store the locations of msigdb category databases
MsigDB.cats = hash(keys=c("Mm.H", "Mm.c2","Mm.c3", "Mm.c4","Mm.c5", "Mm.c6","Mm.c7", "Hs.H", "Hs.c1","Hs.c2","Hs.c3", "Hs.c4","Hs.c5", "Hs.c6","Hs.c7"), 
	values=c("mouse_H_v5p1", "mouse_c2_v5p1", "mouse_c3_v5p1", "mouse_c4_v5p1", "mouse_c5_v5p1", "mouse_c6_v5p1", "mouse_c7_v5p1", 
		"human_H_v5p1",  "human_c1_v5p1","human_c2_v5p1", "human_c3_v5p1", "human_c4_v5p1", "human_c5_v5p1", "human_c6_v5p1", "human_c7_v5p1"))

MsigDB.names = hash(keys = c("H", "c1", "c2", "c3", "c4", "c5", "c6", "c7"), 
	values=c("Hallmark gene sets", "printositional gene sets", "curated gene sets", 
		"motif gene sets", "computational gene sets", "GO gene sets", "oncogenic signatures", "immunologic signatures"))

#GO Analysis object:
 GoAnalysis <- setClass("GoAnalysis", slots = c(go="data.frame", 
 											#dge_list="DGEList", 
 											kegg="data.frame", genes="data.frame"))

go.db = NULL
#source("http://bioconductor.org/biocLite.R")
#biocLite("org.Mm.eg.db")

to.ensembl <- function(genes, db)
{
	ensembls = AnnotationDbi::select(db, keys=names(genes), columns="ENSEMBL", keytype='ALIAS')
	n_no_id = sum(is.na(ensembls$ENSEMBL))
	warn(sprintf("No ENSEMBL id found for %i genes " , n_no_id))
	ensembls = na.omit(ensembls)
	info("Adding diff. expressed column..")
	ensembls["Diff.Expressed"] = genes[match(ensembls$ALIAS,names(genes))]
	n_ambiguous = sum(duplicated(ensembls$ALIAS))
	info(sprintf("More than one ENSEMBL id found for %i genes", n_ambiguous))
	n_ambiguous_other_way = sum(duplicated(ensembls$ENSEMBL))
	info(sprintf("Same ENSEMBL id given to for %i genes.. ", n_ambiguous_other_way))
	warn("Collapsing the ambiguous IDs, keeping only the first ENSEMBL found")
	ensembls = ensembls[!duplicated(ensembls$ALIAS),]
	ensembls = ensembls[!duplicated(ensembls$ENSEMBL),]
	info(sprintf("Found ENSEMBL IDs for %s of %s given genes [%f%%] ", nrow(ensembls), length(genes), nrow(ensembls)/length(genes) * 100))
	return(ensembls)
}

convert.msig <- function(cat, e2s, n.cores=8, output="symbol"){
	info(sprintf("Converting to geneSymbol (multiprocessing using %s cores)", n.cores))
	f = function(x){sort(e2s$symbol[match(x, e2s$gene_id)])}
	# if(output == "ENSEMBL"){
		
	# }else{
	# 	f = function(x){to.ensembl(sort(e2s$symbol[match(x, e2s$gene_id)]))}
	# }
	
	if(n.cores>1){
		mclapply(cat, f, mc.cores=n.cores)
	}else{
		lapply(cat, f )
	}
}

### grab the msigdb gene sets. downloaded from http://bioinf.wehi.edu.au/software/MSigDB/ on sep11
get.msigdb <- function(categories=c("H", "c2", "c5", "c6", "c7"), org="mm10", sort=T)
{
	if(org == "mm10"){
		entrez.to.symbol = toTable(org.Mm.egSYMBOL)
		categories = paste("Mm", categories, sep=".")
	}else{
		if(org == "hg19"){
			library(org.Hs.eg.db)
			categories = paste("Hs", categories, sep=".")
			entrez.to.symbol = toTable(org.Hs.egSYMBOL)
		}else{
			stop(sprintf("Unknown org! --> %s. must be hg19 or mm10", org))
		}
	}

	cat.list = list()
	for(cat in categories){
		file.path = paste0("~/ref/gene_lists/msigDB/", MsigDB.cats[[cat]], ".rds")
		info(sprintf("Loading %s from %s", cat, file.path))
		msig.cat = convert.msig(readRDS(file.path), entrez.to.symbol)
		### keep track of which category each gene set is in:
		names(msig.cat ) = paste(cat, names(msig.cat ), sep=".")
		cat.list = c(msig.cat, cat.list)
		#msig.cat = convert.msig()
	}
	info(sprintf("Loaded %s genesets in %s categories [%s]", length(cat.list), length(categories), org))
	return(cat.list)
}


### wrote this quickly to use for stem-mhc paper, needs to be wrapped in a more general gene_set enrichment function. although tbh it subsumes go analysis.
run.msigdb <- function(genes, background=NULL, 
						org = "mm10",
						working.dir=".",
						analysis.name="MsigDB_enrichment", 
						show_top_terms=30, 
						max.geneset.length=0,
						use.categories=c("H", "c2", "c5", "c6", "c7"), 
						write.excel=T,
						ret.data=F)
{
	library(goseq)
	init.length = length(genes)
	genes = genes[!is.na(genes)]
	genes = unique(genes)
	if(length(genes) < init.length)
	{
		warn(sprintf("Given list shrank from %s to %s genes after removing NAs and duplicates", init.length, length(genes)))
	}
	genes = to.goseq.vector(genes, background)
	
	## something broke this:
	# info("Getting gene lengths..")
	# gene_lengths= getlength(names(genes), org, "geneSymbol")

	gene_lengths = get_gene_lengths(names(genes))

	info("DE genes:")
	de.genes = names(genes)[genes==1]
	#print(de.genes)

	info("Building probability weighting function (PWF)..")
	pwf=nullp(genes, bias.data=gene_lengths)
	
	info("Retrieving mapping from gene names to MSigDB genesets")
	map = get.msigdb(categories=use.categories, org=org)
	if(max.geneset.length > 0){
		init.n = length(map)
		map = map[lapply(map, length) < max.geneset.length]
		warn(sprintf("Removed %s genesets that contain more than %s genes", init.n-length(map), max.geneset.length))
		analysis.name = paste(analysis.name, "maxgenesetlen", max.geneset.length, sep="_")
	}


	info("Running enrichment analysis")
	hits=goseq(pwf, org, "geneSymbol", gene2cat=map)
	
	hits$over_represented_pvalue[hits$over_represented_pvalue==0] <- .Machine$double.eps
	info(sprintf("FDR correcting %s p-values", nrow(hits)))
	hits$over_represented_Qvalue = p.adjust(hits$over_represented_pvalue, method="BH")
	hits = hits[,c(1:2, ncol(hits), 4:ncol(hits)-1)]

	### find the DE genes in each category.
	info("Finding the DE genes in each category")
	
	get_de_genes = function(genes_in_cat){
		e = unlist(genes_in_cat)
		e[toupper(e) %in% toupper(de.genes)]
	}
	
	DEgenes = lapply(map, get_de_genes)
	more_cols = data.frame(unlist(lapply(DEgenes, paste, collapse=", ")), names(DEgenes))
	colnames(more_cols) = c("DE.Genes", "Category_name_check")

	#put go categories in the same order as they come back from goseq() call, ranked by significance.
	more_cols = more_cols[match(hits$category, more_cols$Category_name_check),]
	hits$DE.Genes = more_cols$DE.Genes
	### clean up columns a bit
	hits$under_represented_pvalue <- NULL
	extract.field=function(string,field=1,delim=".", fixed=T) return(strsplit(string,delim, fixed=fixed)[[1]][field])
	hits$MsigDB.Group = unlist(lapply(hits$category, extract.field, 2))
	hits$Org = unlist(lapply(hits$category, extract.field, 1))
	hits$category = unlist(lapply(hits$category, extract.field, 3))

	hits["minus_log10_pval"] = log10(unlist(hits$over_represented_pvalue)) * -1
	hits["minus_log10_Qval"] = log10(unlist(hits$over_represented_Qvalue)) * -1

	i = 1
	plots= list()
	for(cat in sort(unique(hits$MsigDB.Group)))
	{
		info(sprintf("Drawing plot for %s", cat))
		g = hits[as.character(hits$MsigDB.Group)==cat,]
		g = g[1:show_top_terms,]
		if(nrow(g) > 1)
		{
			max_val = max(g$minus_log10_Qval)
			#print(max_val)
			g$category <- factor(g$category, levels=unique(g$category)[order(g$minus_log10_Qval)])
			p = ggplot(g, aes(x=category, y=minus_log10_Qval, fill=minus_log10_pval)) + ylab("Significance (-Log10 Q)") + xlab("") +
			geom_bar(stat='identity', width=0.5)  + coord_flip()  + theme_bw() + 
			theme(axis.text.y = element_text(size=8)) + scale_y_continuous(breaks=seq(0, max_val+1, by=1)) +
			ggtitle(sprintf("%s\nTop %s enriched [GoSeq] gene sets.\nMSigDB group -- %s [%s]", analysis.name, show_top_terms, MsigDB.names[[cat]], cat)) +
			scale_fill_gradientn( colours=colorRampPalette(brewer.pal(9,"Blues"))(50), guide=FALSE) +
			theme(plot.title = element_text(lineheight=.8, face="bold", size=10)) + facet_wrap(~MsigDB.Group)
			#print(p)
		}else{
			warn(sprintf("No enriched %s terms!", category))
		}
		plots[[i]] <- p
		i = i + 1
	}
	pdf.dest = paste0(working.dir, "/", sprintf("%s.pdf", analysis.name))
	info(sprintf("Rendering to %s", pdf.dest))
	pdf(pdf.dest, width=12, height=8)
	invisible(lapply(plots, print))
	dev.off()

	table.dest = paste0(working.dir, "/", sprintf("%s.txt", analysis.name))
	info(sprintf("Saving results to %s", table.dest))
	write.table(hits, file=table.dest, sep="\t", quote=F, row.names=F)
	
	#### write a formatted table to excel spreadsheet.
	if(write.excel){
		excel.dest = paste0(working.dir, "/", sprintf("%s.xlsx", analysis.name))
		library(xlsx)
		options(java.parameters = "-Xmx16000m")
		if (file.exists(excel.dest)) file.remove(excel.dest)
		wb <-createWorkbook(type="xlsx")
		TITLE_STYLE     <- CellStyle(wb) + Font(wb, heightInPoints=16, isBold=TRUE, underline=1, color = "blue")
		SUB_TITLE_STYLE <- CellStyle(wb) + Font(wb, heightInPoints=14, isItalic=TRUE, isBold=FALSE, color = "blue")

		TABLE_ROWNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) + Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") + 
								Border(color="black", position=c("RIGHT"), 
								           pen=c("BORDER_THIN")) 
		TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) +
		    Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
		    Border(color="black", position=c("TOP", "BOTTOM"), 
		           pen=c("BORDER_THIN", "BORDER_THICK")) 
	    
        CELL_STYLE <-  CellStyle(wb) + Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER")

		xlsx.addTitle<-function(sheet, rowIndex, title, titleStyle){
				  rows <-createRow(sheet,rowIndex=rowIndex)
				  sheetTitle <-createCell(rows, colIndex=1)
				  setCellValue(sheetTitle[[1,1]], title)
				  setCellStyle(sheetTitle[[1,1]], titleStyle)
			}


		for(cat in unique(hits$MsigDB.Group))
		{
			sheetname = MsigDB.names[[cat]]
			sheet <- createSheet(wb, sheetName = sheetname)
			xlsx.addTitle(sheet, rowIndex=2, 
			              title=analysis.name,
			              titleStyle = TITLE_STYLE)
			xlsx.addTitle(sheet, rowIndex=4, 
				      title=sprintf("Enriched gene sets [GoSeq] using MSigDB -- %s [%s]", MsigDB.names[[cat]], cat),
				      titleStyle = SUB_TITLE_STYLE)
			df  = hits[hits$MsigDB.Group==cat,]
			rownames(df) = make.names(df$category, unique=T)
			df$category <- NULL
			addDataFrame(df, sheet, startRow=6, startColumn=1, 
					 colStyle = CELL_STYLE,
		             colnamesStyle = TABLE_COLNAMES_STYLE,
		             rownamesStyle = TABLE_ROWNAMES_STYLE)
			setColumnWidth(sheet, colIndex=c(1, 6), colWidth=60)

			#write.xlsx2(, file=, sheetName=sheetname, append=TRUE,row.names=F)
		}
		saveWorkbook(wb, excel.dest)
		jgc()
	}
	
	if(ret.data){
		return(hits)
	}
}


comparative.go.plot <- function(go1.results, go2.results, names=c("GO_1", "GO_2"), 
	use.column=NULL, ## if this is set, only plot terms with a 'T' in the column with this name
	use.top.terms=0, 
	fdr.cutoff=0.1, 
	pdf.output="Comparative_GO.pdf",
	width=12,
	height=9,
	#down.cols=brewer.pal(9, "Blues"),
	#up.cols=brewer.pal(9, "Reds"), 
	cols=rev(brewer.pal(9, "RdBu")),
	breaks.by=1,
	hjust.down=2,
	hjust.up=-0.1,
	text.size.term=18,
	extra.space.up = 3,
	extra.space.down = 3)
{

	treatment_go = read.delim(go1.results)
	ctl_go = read.delim(go2.results)

	ctl_go = ctl_go[ctl_go$ontology=="Biological.Process",]
	treatment_go = treatment_go[treatment_go$ontology=="Biological.Process",]
	
	if(is.null(use.column)){
		if(fdr.cutoff > 0){
			ctl_go = ctl_go[ctl_go$over_represented_Qvalue < fdr.cutoff,]
			treatment_go = treatment_go[treatment_go$over_represented_Qvalue < fdr.cutoff,]
		}else{
			if(use.top.terms > 0)
			{
				ctl_go = ctl_go[1:use.top.terms,]
				treatment_go = treatment_go[1:use.top.terms,]
			}else{
				stop("Provide an FDR cutoff > 0 or set use.top.terms > 0")
			}
		}
	}else{
		info(sprintf("Using only terms with a \'T\' in \'%s\' column", use.column))
		print(head(ctl_go))
		print(table(ctl_go$use))
		print(table(treatment_go$use))
		ctl_go = ctl_go[!is.na(ctl_go$use),]
		ctl_go = ctl_go[which(ctl_go$use), ]
		
		treatment_go = treatment_go[!is.na(treatment_go$use),]
		treatment_go = treatment_go[which(treatment_go$use), ]

		#print(ctl_go[, c("term", "over_represented_Qvalue")])
	}
	
	
	n.terms.up = nrow(ctl_go)
	n.terms.down = nrow(treatment_go)

	info(sprintf("%s terms UP, %s DOWN", n.terms.up, n.terms.down))

	x = rbind(ctl_go, treatment_go)

	x["minus_log10_Qval"] = log10(as.numeric(x$over_represented_Qvalue)) * -1
	max_val = ceiling(max(x$minus_log10_Qval))
	x$Condition = c(rep(names[1], n.terms.up), rep(names[2], n.terms.down))

	down.terms = (n.terms.up+1):(n.terms.up + n.terms.down)
	x$minus_log10_Qval[down.terms] = -1 * x$minus_log10_Qval[down.terms]

	brk.size = 1 + (max_val %/% 10)
	ggsave(ggplot(x, aes(x=reorder(GO_TERM_check, abs(minus_log10_Qval)), y=minus_log10_Qval, fill=minus_log10_Qval, label=term)) + 
		ylab("Significance: -Log10(Q)") + xlab("") +
  		geom_bar(stat='identity', width=0.5)  + coord_flip()  + theme_bw() + 
  		scale_y_continuous(breaks=seq(-max_val-extra.space.down, max_val+extra.space.up, by=breaks.by), limits = c(-max_val-extra.space.down, max_val+extra.space.up), expand=c(0,0)) +
  		theme(plot.title = element_text(lineheight=.8, face="bold", size=10)) + geom_text(data=subset(x, minus_log10_Qval>0), size=text.size.term, hjust=hjust.up) + 
  		geom_text(data=subset(x, minus_log10_Qval<0),  size=text.size.term, hjust=hjust.down) + 
  		theme(axis.line.x = element_line(color="gray15", size = 0.75, lineend="round"), 
  			axis.text.x = element_text(size=12), panel.grid.major = element_blank(), 
  			panel.grid.minor = element_blank(), panel.border= element_blank(), 
  			panel.background =  element_blank(), axis.ticks.y=element_blank(), 
  			axis.text.y=element_blank()) + 
  	scale_fill_gradientn("Significance", colors = cols, guide = FALSE), filename=pdf.output, width=width, height=height)

}


# draw the same plot as in 'marker.enrichment.test'
# except combine the up and downregulated gene sets 
combined.plot <- function(up, down, filename, max.value=10)
{	
	info(sprintf("Saving combined plot to %s", filename))
	x = rbind(up, down)
	x$type =rep(c("Up", "Down"), each=nrow(up))
	x$logp = -log10(x$p)
	x$logp[x$type=="Down"] <- x$logp[x$type=="Down"] * -1
	
	if(!all(is.finite(x$p))){
		warn("p val column contains inifinite values (NaNs)! skipping plot")
		return()
	}
	
	if(sum(x$p)==nrow(x))
	{
		warn("No significant hits, not drawing plot")
		return()
	}
	print(x)
	g=ggplot(x, aes(x=MarkerGeneSet, y=logp, fill=type)) + geom_bar(stat="identity", position="dodge", width=0.5) + theme_bw() + 
		xlab("") + geom_hline(yintercept=-log10(0.05), linetype="dotted", color="firebrick2") + 
		geom_hline(yintercept=0,  color="black", size=1) + 
		geom_hline(yintercept=log10(0.05), linetype="dotted", color="firebrick2") + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Prob. of up or down regulation (-log p)") + 
		scale_fill_manual("", values=default.cols(2)) + coord_cartesian(ylim=c(-max.value,max.value))
	ggsave(g, filename=filename, width=9, height=7)
}





convert_entrez_to_gene_symbols <- function(entrez.ids, entrez.to.symbol)
{
	entrez.rows = match(entrez.ids, unlist(entrez.to.symbol$gene_id))
	#cat(sprintf("Found entrez ids for %s genes \n", length(entrez.rows)))
	tb = entrez.to.symbol[entrez.rows,]
	#print(tb)
	return (unlist(tb$symbol))
}

# a utility to read in a table of markers and write an updated table
# with the GO annotations:
annotate.markers.file <- function(filename)
{
	df = read.delim(filename)

	df.ann = annotate.markers(df)
	#df.ann["GENE_SYMBOL"] = rownames(df.ann)
	df.ann = cbind(GENE_SYMBOL=rownames(df.ann), df.ann)
	write.table(df.ann, 
			file="with_annotated_markers.txt", 
			sep="\t", 
			col.names=T,
			row.names=F)
}

# given a dataframe with gene names as the rownames
# add two columns, one with important go terms and
# one with all go terms.
# parallel gives weird results, maybe the GO thing isn't parallelisable.
annotate.markers <- function(data) #, parallel=T)
{
	print("Received data frame:")
	print(dim(data))
	# if(parallel)
	# {
	# 	terms = unlist(mclapply(rownames(data), function(x) get.gene.annotations(x, all.gene.names=rownames(data), verbose=T), mc.cores=detectCores()))
	# }else{
		terms = sapply(rownames(data), function(x) get.gene.annotations(x, all.gene.names=rownames(data), verbose=T))
	# }
	print("Calculated go terms:")
	print(terms)
	print("Size:")
	print(dim(terms))
	data["Important.Go.Terms"] = unlist(terms[1,])
	data["Go.Terms"] = unlist(terms[2,])
	print("Results data frame:")
	print(dim(data))
	return (data)
}


get.alias <- function(gene.name)
{
	library(GO.db)
	library(org.Mm.eg.db)
	library(DBI)
	dbCon <- org.Mm.eg_dbconn()
	# write your SQL query
	sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
	# execute the query on the database
	aliasSymbol <- dbGetQuery(dbCon, sqlQuery)
	# subset to get your results
	aliasSymbol[which(aliasSymbol[,2] == gene.name),5]
}


# for a given gene, return GO annotations.
# useful info: http://www.bioconductor.org/help/workflows/annotation-data/
# elements of returned value:
# 1. terms matching the important.terms argument concatenated into a string
# 2. all terms concatenated into a string
# 3. all terms as a list

## MOUSE ONLY 

get.gene.annotations <- function(gene.name, 
									id.type="GENE_SYMBOL", 
									all.gene.names=NULL, # keep track of index
									verbose=F,
									important.terms = c("cell surface", "regulation of transcription", "membrane"))
{
	
	if(verbose)
	{
		if(!is.null(all.gene.names))
		{
			total.genes = length(all.gene.names)
			current.gene = which(all.gene.names==gene.name)
			cat(sprintf("\nLooking for info for gene: %s 		[%i of %i]\n", gene.name, current.gene, total.genes))
		}else{
			cat(sprintf("\nLooking for info for gene: %s \n", gene.name))

		}
	}

	library(GO.db)
	library(org.Mm.eg.db)

    rval <- tryCatch(
        {
		    #print(gene.name)
			# cat(sprintf("Converting to ENTREZ..\n"))
			# entrez = unlist(gene.sym.to.entrez[gene.name])
			# #annots <- select(org.Mm.eg.db,keys=entrez,columns="SYMBOL",keytype="ENTREZID")
			#cat(sprintf("Getting GO-to ENTREZ..\n"))
			
			# use sql to get alias table and gene_info table (contains the symbols)
			# first open the database connection
			gene.name = get.alias(gene.name) ### make sure we have the right symbol.
			go.id = AnnotationDbi::select(org.Mm.eg.db, gene.name, "GO", "SYMBOL")
			


			#print(go.id[entrez[[1]]]$Ontology)

			go.terms = AnnotationDbi::select(GO.db, go.id$GO, "TERM", "GOID")
			important = go.terms[grep(paste(important.terms, collapse="|"),go.terms$TERM), ]
			
			if(verbose)
			{
				cat(sprintf("Returning %i GO terms, %i of match %s .. \n\n", length(unlist(go.terms$TERM)), length(unlist(important$TERM)), paste(important.terms, collapse=" OR ")))
			}

			as.long.string = paste(unlist(go.terms$TERM), sep=";", collapse="; ")
			important.as.long.string = paste(unlist(important$TERM), sep=";", collapse="; ")
			rval =  (c(important.as.long.string, as.long.string, go.terms))
        },
        error=function(cond) {
            print(sprintf("ERROR: Gene lookup failed! [%s]", gene.name))
            print(cond)
            info.string = "GO Lookup Failed!"
            x = data.frame(cbind(info.string, info.string))
            colnames(x) = c("GOID", "TERM")
            rval =  c(c(NA),c(NA),x)
            
        },
        warning=function(cond) {
            #ignore warnings
        },
        finally={
        	#don't do anything
        }
    )    
	#print("Returning:")
	#print(rval)
	#print("Size:")
	#print(length(rval))
    return(rval)
}

### convert genes to the format GOSeq expects:
### In order to perform a GO analysis of your RNA-seq data, goseq only requires a simple named
### vector, which contains two pieces of information.
### 1. Measured genes: all genes for which RNA-seq data was gathered for your experiment. Each
### element of your vector should be named by a unique gene identifier.
### 2. Differentially expressed genes: each element of your vector should be either a 1 or a
### 0, where 1 indicates that the gene is differentially expressed and 0 that it is not.
to.goseq.vector <- function (gene_list, background=NULL){
	if(is.null(background)){
		warn("Using some random gene list as background. Don't use this for final analysis!")
		if(!check_yes_no()){
			stop()
		}
		background = as.character(read.delim(
			"/Users/ahaber/Dropbox/Postdoc/Projects/GutCircuits/Analysis/Datasets/15allyear_GutAtlas_SS2/May_merge/Atlas_and_MHC_paper/Tuft_Cells_Figure_4/GO_analysis/background_genes.txt"
			)$x)
	}
	background = background[!toupper(background) %in% toupper(gene_list)]
	goseq.vector = c(rep(1, length(gene_list)), rep(0, length(background)))
	names(goseq.vector) = c(gene_list, background)
	print(table(goseq.vector))
	print(head(goseq.vector))
	return(goseq.vector)
}

### main wrapper function to call GO. just converts the genes to the right format (vector) and makes a directory to work in.
go.analysis <- function(gene_list, background=NULL, analysis.name="GO.Analysis", annotation="mm10")
{
	init.dir = getwd()
	dir.create(analysis.name)
	setwd(analysis.name)
	goseq.vector = to.goseq.vector(gene_list, background)
    go = run_go(goseq.vector, analysis_group_name=analysis.name, show_top_terms=25, annotation=annotation)
    setwd(init.dir)
}

# the input to the function is assumed to be an edger DGEList
# object, that specifies which genes are differentially expressed
# in the population of interest.
edger_to_go_vector <- function(markers_table, p = 0.05)
{
	cat(sprintf("Input object: \n"))
	print(head(markers_table))
	genes=as.integer(p.adjust(markers_table$table$PValue[edger_tested$table$logFC!=0], method="BH") < p)
	n_genes = as.matrix(table(genes))[2]
	cat(sprintf("Using %i genes upregulated(p < %f) for GO enrichment analysis.. [Total genes=%i] \n", n_genes, p, length(genes)))
	gene_symbols = row.names(markers_table$table[edger_tested$table$logFC!=0,])
	names(genes) = gene_symbols
	return (genes)
}

scde_to_go_vector <- function(markers_table, p =0.05, use.top=0)
{
	# cat(sprintf("Input object: \n"))
	# print(head(markers_table))
	if(use.top > 0)
	{
		warn(sprintf("Using top %s markers (ignoring p-value threshold)", use.top))
		genes=as.integer(1:nrow(markers_table) < use.top)
	}else{
		if("Mean.Min.LogFC" %in% colnames(markers_table))
		{
			info("Using genes DE against both background and similar cluster!")
			genes=as.integer(markers_table$Lower.bound.background > 1 & 
				markers_table$p.background < p & markers_table$p.sim.clust)
		}else
		{
			genes=as.integer(markers_table$Lower.bound > 1 & markers_table$p < p)
		}
		
		
		
	}
	n_genes = as.matrix(table(genes))[2]
	cat(sprintf("Using %i genes upregulated(p < %f) for GO enrichment analysis.. [Total genes=%i] \n", n_genes, p, length(genes)))
	names(genes) = markers_table$GENE_SYMBOL
	genes = genes[sort(names(genes))]
	
	return (genes)
}

# called by the single cell pipeline to run GO analysis
# on each cluster of cells using the SCDE markers.
run.go.analysis <- function(counts, clusters, 
	p.threshold=0.05, 
	use.top.genes=0)
{

	info(sprintf("Running GO analysis on %i clusters..\n", length(unique(clusters))))
	info(sprintf("Use top genes     -> %s ", use.top.genes))
	info(sprintf("p-value threshold -> %s ", p.threshold))
	cluster.names = unique(clusters)
	initial_dir = getwd()

	
	for(i in 1:length(cluster.names))
    {
        current.cluster = as.character(cluster.names[i])
        if(use.top.genes > 0)
        {
        	all.genes.file = paste(current.cluster, paste("markers_", current.cluster, ".txt", sep=""), sep="/")
    	}else{
    		all.genes.file = paste(current.cluster, paste("all_genes_", current.cluster, ".txt", sep=""), sep="/")
    	}
        
        info(sprintf("Reading DE genes for cluster %i of %i from %s..", i, length(cluster.names), all.genes.file))

        if (isTRUE(file.exists(all.genes.file))) 
        {
        	markers = read.delim(all.genes.file)
        	setwd(current.cluster)
	        goseq.vector = scde_to_go_vector(markers, p=p.threshold, use.top=use.top.genes)
	        gl = NULL # should only load gene lengths once
	        go = run_go(goseq.vector, analysis_group_name=current.cluster, show_top_terms=25, annotation="mm10")
	        in.current = clusters == current.cluster
	        if(!is.null(go))
	        {
        		#visualise.pathways(go, counts, cluster.ids=in.current)
        	}else{
        		warn("No go analysis, cannot visalise KEGG pathways..")
        	}
        	setwd(initial_dir)
        }else{
        	warn(sprintf("No genes found for %s!", current.cluster))
        }
    }
    info("GO analysis done!")
}

test.go <- function(scde.markers, p.threshold=0, use.top.genes=20)
{
	markers = read.delim(scde.markers)
    goseq.vector = scde_to_go_vector(markers, p=p.threshold, use.top=use.top.genes)
    print(head(goseq.vector))
    go = run_go(goseq.vector, analysis_group_name="TEST", show_top_terms=25)
}

# functions to run gene ontology (GO) analysis a list of (usually DE) genes
run_go = function(genes, 
				  counts = NULL,
				  annotation="mm10", 
				  gene_id="ensGene", 
				  analysis_group_name="Unknown", 
				  go.cat.only="BP",						#can be "CC", "BP", or "MF" for ; Cellular Components, Biological
														#Processes and Molecular Functions. or any combination of them.
				  produce_plot=TRUE,
				  plot_destination=NULL,
				  show.q=0.5,
				  min.genes=5,
				  show_top_terms=50,
				  run_kegg=T,
				  run_msigdb=T,
				  lookup.names=T)
{
	info("Loading GO analysis libraries")
	library(goseq)
	library(plyr)
	library(biomaRt)
	library(GO.db)
	## see bottom of pg 13 of goseq manual: https://bioconductor.org/packages/devel/bioc/vignettes/goseq/inst/doc/goseq.pdf
	go.cat.only = paste0("GO:", go.cat.only)
	info(sprintf("Beginning GO analysis of %s [GoSeq] ", analysis_group_name))
	info(sprintf("Annotation: %s ",annotation))
	info(sprintf("Gene Id: %s ",gene_id))


	if(annotation=="mm10"){
		# Mouse
		library(org.Mm.eg.db)
		library("Mus.musculus")
		entrez.to.symbol = toTable(org.Mm.egSYMBOL)
		entrez.to.go <- as.list(org.Mm.egGO)
		gene.sym.to.entrez <- as.list(org.Mm.egALIAS2EG)
		# for kegg pathways
		org.id = "mmu"
		db = Mus.musculus
	}else{
		if(annotation=="hg19")
		{
			# Human
			library(Homo.sapiens)
			library(org.Hs.eg.db)
			entrez.to.symbol = toTable(org.Hs.egSYMBOL)
			entrez.to.go <- as.list(org.Hs.egGO)
			gene.sym.to.entrez <- as.list(org.Hs.egALIAS2EG)
			# for kegg pathways
			org.id = "hsa"
			db = Homo.sapiens
		}else{
			stop(sprintf("Unknown annotation! --> %s", annotation))
		}
	}
	
	id = "ENSEMBL" #ENSEMBLTRNS
	cat(sprintf("Converting to %s", id))
	
	# if(!is.null(counts))
	# {
	# 	# use refseq ids in the transcript_id(s) column to convert to ENSEMBLTRANS
	# 	cat(sprintf(" using RefSeq ids.. \n"))
	# 	f = function(x){if(grepl(",",x)){strsplit(as.character(x),",")[[1]][1]}else{as.character(x)}}
	# 	refseq_ids= unlist(lapply(counts$transcript_id.s., f))
	# 	ensembls = AnnotationDbi::select(db, keys=refseq_ids, columns=id,keytype='REFSEQ')
	# 	# watch out for duplicates:
	# 	# collapse the ambiguous values so order is (i hope) preserved
	# 	cat(sprintf("Collapsing duplicates and counting missing IDs.. \n"))
	# 	ensembls = ddply(ensembls, "REFSEQ", summarize, ENSEMBL = paste(ENSEMBL, collapse="&"))
		
	# }else{
		# use genesymbol (rownames)
	info("using GENE_SYMBOL")
	#}

	ensembls = to.ensembl(genes, db)

	if(nrow(ensembls) < min.genes)
	{
		error("Not enough genes found to run analysis!")
		return (NULL)
	}

	genes = ensembls$Diff.Expressed
	names(genes) = ensembls$ENSEMBL

	info("Getting gene lengths..")
	gene_lengths= getlength(names(genes), annotation, gene_id)

	info("Building probability weighting function (PWF)..")
	pwf=nullp(genes,bias.data=gene_lengths)
	
	# PWF - probability weighting function
	# http://www.bioconductor.org/packages/release/bioc/vignettes/goseq/inst/doc/goseq.pdf 
	# info("Compensating for gene length using PWF:")
	# print(head(pwf))

	info("Retrieving mapping from gene names to GO terms")
	go_map=getgo(names(genes), annotation, gene_id)

	info(sprintf("Running GO analysis [cats=%s]", paste(go.cat.only, collapse=", ")))
	go=goseq(pwf, annotation, gene_id, gene2cat=go_map, test.cats=go.cat.only)
	
	print(head(go))
    cats = hash(keys=c("BP", "CC","MF"), values=c("Biological.Process", "Cellular.Component", "Molecular.Function"))
	types.used = cats[[go.cat.only]]
	info(sprintf("GO categories used: %s", paste(types.used, collapse=", ")))
	go$ontology = mapvalues(go$ontology, hash::keys(cats), hash::values(cats))
	
	info("Generating reverse mapping (GO -> genes)")
	go2genes=goseq:::reversemapping(go_map)
	info("Finding the DE genes in each category")
	
	convert_to_de_gene_symbols = function(ens_ids){
		e = unlist(ens_ids)
		symbols = ensembls[ensembls$ENSEMBL %in% ens_ids, ]
		symbols = symbols[symbols$Diff.Expressed==1,]
		return (symbols$ALIAS)
	}
	
	DEgenes = lapply(go2genes, convert_to_de_gene_symbols)
	more_cols = data.frame(unlist(lapply(DEgenes, paste, collapse=", ")), names(DEgenes))
	colnames(more_cols) = c("Enriched.Genes", "GO_TERM_check")

	#put go categories in the same order as they come back from goseq() call, ranked by significance.
	more_cols = more_cols[match(go$category, more_cols$GO_TERM_check),]
	go$Enriched.Genes = more_cols$Enriched.Genes
	go$GO_TERM_check = more_cols$GO_TERM_check

	# print(head(go, n=12))
	# #let's annotate the GO categories
	# cats <- go$category
	# terms <- stack(lapply(mget(cats, GOTERM, ifnotfound=NA), Term))
	# go$Term <- with(go, terms$values[match(terms$ind, go$category)] )
	# allGos <- stack(getgo(rownames(topGenes$table), 'bosTau4', 'ensGene')) # so here I pull the GO terms for every gene that is regulated.

	# if(!is.null(go.cat.only))
	# {
	# 	info("Restricting hits to those in categories:")
	# 	info(paste(go.cat.only, collapse=", "))
	# 	n.init = nrow(go)
	# 	go = go[go$ontology %in% go.cat.only,]
	# 	info(sprintf("Reduced from %s to %s GO terms", n.init, nrow(go)))
	# 	#print(table(go$ontology))
	# }

	### Do multiple hypothesis testing.. method="BH" is the same as method="fdr"
	### so resulting values are FDR-corrected p values (Q values).
	go[go==0] <- .Machine$double.eps
	go$over_represented_Qvalue = p.adjust(go$over_represented_pvalue, method="BH")
	###

	n.sig = sum(go$over_represented_Qvalue < show.q)
	info(sprintf("Found %s pathways significantly enriched (FDR < %s)", n.sig, show.q))
	

	# draw a plot of enriched go terms:
	if(produce_plot){
		go["minus_log10_pval"] = log10(unlist(go$over_represented_pvalue)) * -1
		go["minus_log10_Qval"] = log10(unlist(go$over_represented_Qvalue)) * -1
		info("Drawing enriched GO terms plot")

		if(is.null(plot_destination)){
			plot_destination = paste(analysis_group_name, "GO_terms.pdf", sep="_")
		}
		pdf(plot_destination, width=12, height=8)
		
		#one plot per go-ontology category, there are 3. BP is the most relevant, do 2 first.
		for(category in types.used)
		{
			info(sprintf("Building GO plot for %s..", category))
			g = head(go[go$ontology == category,], n=show_top_terms)
			if(nrow(g) > 1)
			{
				max_val = max(g$minus_log10_Qval)
				#print(max_val)
				g$term <- factor(g$term, levels=rev(unique(g$term)))
				p = ggplot(g, aes(x=term, y=minus_log10_Qval, fill=minus_log10_pval)) + ylab("Significance (-Log10 Q)") + xlab("") +
				geom_bar(stat='identity', width=0.5)  + coord_flip()  + theme_bw() + 
				theme(axis.text.y = element_text(size=8)) + scale_y_continuous(breaks=seq(0, max_val+1, by=1)) +
				ggtitle(gsub(".", " ", paste("Top", show_top_terms, category,"GO terms in", analysis_group_name, "cluster"), fixed=T)) +
				scale_fill_gradientn( colours=colorRampPalette(brewer.pal(9,"Blues"))(50), guide=FALSE) +
				theme(plot.title = element_text(lineheight=.8, face="bold", size=10))
				print(p)
			}else{
				warn(sprintf("No enriched %s terms!", category))
			}
		}
		dev.off()
	}
	
	# cat(sprintf("Results: \n"))
	go.table=paste(analysis_group_name, "GO_terms.txt", sep="_")
	write.table(go, file=go.table, sep="\t", quote=F, row.names=F)

	DE.gene.symbols = unlist(ensembls[which(ensembls$Diff.Expressed==1), "ALIAS"])
	print(DE.gene.symbols)
	
	###TODO fold the following code for KEGG pathway analysis (which works)
	KEGG = data.frame()

	# if(run_kegg){
	# 	cat(sprintf("Running KEGG pathway analysis.. \n"))
	# 	KEGG=goseq(pwf,annotation, gene_id,test.cats="KEGG")
	# 	if (annotation=="mm10")
	# 	{
	# 		org = "mmu"
	# 	}else if(annotation == "hg19")
	# 	{
	# 		org = "hsa"
	# 	}else
	# 	{
	# 		warn(sprintf("Cannot look up KEGG pathway names for unknown organism: %s" , annotation))
	# 	}
	# 	if(lookup.names)
	# 	{
	# 		use.kegg.db = TRUE # much faster, but could possibly be out of date, see:
	# 		# http://bioconductor.org/packages/release/data/annotation/manuals/KEGG.db/man/KEGG.db.pdf
	# 		names = list()
	# 		genes.entrez = list()
	# 		genes.de = list()
	# 		if(use.kegg.db)
	# 		{
	# 			library(KEGG.db)
	# 			xx <- as.list(KEGGPATHNAME2ID)
	# 			yy <- as.list(KEGGPATHID2EXTID)
	# 			# print(head(yy))

	# 			for(i in 1:nrow(KEGG))
	# 			{
	# 				#info(sprintf("Naming pathway %i of %i..", i, nrow(KEGG)))
	# 				id = KEGG$category[i]
	# 				name = names(xx)[grep(id, xx)]
	# 				id = paste(org.id, id, sep="")
	# 				names[[i]] = name
	# 				#cat(sprintf("ID: %s\n", id))
	# 				entrez.symbols = unlist(yy[id])
	# 				#cat("Genes: \n")
	# 				#print(entrez.symbols)
	# 				gene.symbols = convert_entrez_to_gene_symbols(entrez.symbols, entrez.to.symbol)
	# 				# print("Converting")
	# 				# print(head(entrez.symbols))
	# 				# print("to")
	# 				# print(head(gene.symbols))
	# 				de = which(gene.symbols %in% DE.gene.symbols)
	# 				# cat(sprintf("Found gene symbols for %s de genes\n", length(de)))
	# 				de.string = paste(gene.symbols[de], collapse=",")
	# 				genes.de[[i]] = de.string
	# 				# print("DE:")
	# 				# print( de.string)
	# 			}
	# 		}else{
	# 			library(RCurl)
	# 			for(i in 1:nrow(KEGG))
	# 			{
	# 				info(sprintf("Naming pathway %i of %i..", i, nrow(KEGG)))
	# 				id = KEGG$category[i]
	# 				name = getURL(paste("http://togows.dbcls.jp/entry/pathway/", org, id, "/name", sep=""))
	# 				names[[i]] = name
	# 			}
	# 		}
			
	# 		KEGG["Name"] = unlist(names)
	# 		#KEGG["Genes"] = unlist(genes.entrez)
	# 		KEGG["DE.Genes"] = unlist(genes.de)
	# 		KEGG[KEGG==0] <- .Machine$double.eps
	# 		KEGG$over_represented_Qvalue = p.adjust(KEGG$over_represented_pvalue, method="BH")

	# 		kegg.table=paste(analysis_group_name, "KEGG_pathways.txt", sep="_")
	# 		write.table(KEGG, file=kegg.table, sep="\t", quote=F, row.names=F)

	# 		if(produce_plot)
	# 		{
				
	# 			KEGG = head(KEGG, n=show_top_terms)
	# 			KEGG["minus_log10_pval"] = log10(as.numeric(KEGG$over_represented_pvalue)) * -1
	# 			KEGG["minus_log10_Qval"] = log10(as.numeric(KEGG$over_represented_Qvalue)) * -1

	# 			max_val = max(KEGG$minus_log10_Qval)
	# 			# print(max_val)
	# 			KEGG$Name <- factor(KEGG$Name, levels=rev(unique(KEGG$Name)))
	# 			#print(KEGG)
				
	# 			p = ggplot(KEGG, aes(x=Name, y=minus_log10_Qval, fill=minus_log10_Qval)) + ylab("Significance (-Log10 Q)") + xlab("") +
	# 					geom_bar(stat='identity', width=0.5)  + coord_flip()  + theme_bw() + 
	# 					theme(axis.text.y = element_text(size=8)) + scale_y_continuous(breaks=seq(0, max_val+1, by=1)) +
	# 					ggtitle(gsub(".", " ", paste("Top", show_top_terms, "enriched KEGG pathways in", analysis_group_name, "cluster"), fixed=T)) +
	# 					scale_fill_gradientn( colours=colorRampPalette(brewer.pal(9,"Blues"))(50), guide=FALSE) +
	# 					theme(plot.title = element_text(lineheight=.8, face="bold", size=10))
	# 					plot_destination = paste(analysis_group_name, "KEGG_pathways.pdf", sep="_")
	# 					ggsave(p, file=plot_destination, width=12, height=8 )
	# 		}
	# 	}

	# }else{
	# 	KEGG=NULL
	# }

	
	rval = GoAnalysis(go=go, kegg=KEGG, genes=ensembls)
	return (rval)
}

# use the pathview package to visusalise upregulated pathways
visualise.pathways <- function( go, 
								counts, 
								cluster.ids = NULL, #include a boolean vector for the cells whose gene expression to look at
								n=5, 
								pdf.output=T, 
								pdf.name="kegg.pathways.pdf", 
								output.dir = "kegg.pathway", 
								dataset.name="gut.rnaseq.")
{
	info(sprintf("Visualising top %i upregulated pathways..", n))
	library(pathview)	
	
	gene.data = rowMeans(counts)
	if(!is.null(cluster.ids))
	{
		info("Using cluster ids: ")
		print(table(cluster.ids))
		counts = t(scale(t(counts)))
	}

	if(pdf.output)
	{
		pdf(pdf.name, width=12, height=8)
	}

	initial.dir = getwd()
	dir.create(output.dir, showWarnings = FALSE)
	#pipeline.output.dir = file_path_as_absolute(output.dir)
	setwd(output.dir)
	for(i in 1:n)
	{
		info(sprintf("Pathway %i, %s..", i, go@kegg$Name[i]))
		# graphviz:
		pv.out <- pathview(gene.data = gene.data, pathway.id = go@kegg$category[i], 
			species = "hsa", out.suffix = dataset.name, kegg.native = F, gene.idtype="GENENAME")
		# kegg:
		pv.out <- pathview(gene.data = gene.data, pathway.id = go@kegg$category[i], 
			species = "hsa", out.suffix = dataset.name, gene.idtype="GENENAME")
	}
	if(pdf.output)
	{
		dev.off()
	}
	setwd(initial.dir)
	info("Done.")
}


### added this in march 2018 after the getlengths function no longer worked 
get_gene_lengths <- function(geneSymbols)
{

	info("Loading mm10 libraries")
    library(TxDb.Mmusculus.UCSC.mm10.knownGene) 
    library(org.Mm.eg.db)
    txsByGene = transcriptsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by='gene')    
    symbol_to_ucsc_id    = unlist(  mget(geneSymbols[ geneSymbols %in% keys(org.Mm.egSYMBOL2EG) ],org.Mm.egSYMBOL2EG) )
    lengthData=median(width(txsByGene)) ## from the goseq manual, section 5.3 http://www.bioconductor.org/packages/3.7/bioc/vignettes/goseq/inst/doc/goseq.pdf
   
    info("Getting gene lengths (takes ~20 seconds)")
    lengths = rep(NA, length(geneSymbols))
    names(lengths) = geneSymbols
    for(i in 1:length(geneSymbols)){
    	ucsc_id = symbol_to_ucsc_id[geneSymbols[i]]
    	lengths[i] <- lengthData[ucsc_id]
    }
    info("Done!")
    return(lengths)
}


# get the lengths of the genes from the GTF and Fasta files used
# to estimate expression abundances:

# # NOT USED OR TESTED
# get_lengths = function( GTFfile = "something.gtf", 
# 						FASTAfile = "something.fa"){
# 	library(GenomicRanges)
# 	library(rtracklayer)
# 	library(Rsamtools)


# 	#Load the annotation and reduce it
# 	GTF <- import.gff(GTFfile, format="gtf", genome="mm10", asRangedData=F, feature.type="exon")
# 	grl <- reduce(split(GTF, elementMetadata(GTF)$gene_name))
# 	reducedGTF <- unlist(grl, use.names=T)
# 	elementMetadata(reducedGTF)$gene_name <- rep(names(grl), elementLengths(grl))

# 	#Open the fasta file
# 	FASTA <- FaFile(FASTAfile)
# 	open(FASTA)

# 	#Add the GC numbers
# 	elementMetadata(reducedGTF)$nGCs <- letterFrequency(getSeq(FASTA, reducedGTF), "GC")[,1]
# 	elementMetadata(reducedGTF)$widths <- width(reducedGTF)

# 	#Create a list of the ensembl_id/GC/length
# 	calc_GC_length <- function(x) {
# 	    nGCs = sum(elementMetadata(x)$nGCs)
# 	    width = sum(elementMetadata(x)$widths)
# 	    c(width, nGCs/width)
# 	}
# 	output <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_name), calc_GC_length))
# 	colnames(output) <- c("Length", "GC")

# 	#write.table(output, file="GC_lengths.tsv", sep="\t")
# 	return (output$Length)
# }



