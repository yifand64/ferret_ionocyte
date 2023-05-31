
### given a table of SS2 genes, convert any genes for which there is 
### mismatch between SS2/10x names, and we can replace the SS2 name with the
### 10x one. Mismatches where there is no replacement (because there's no ensembl)
### are left alone.
convert.table.ss2.to.10x <- function(tab.ss2)
{
	tab = read.delim("~/ref/tables/10x_SS2_genes.txt")
	for(sig in colnames(tab.ss2))
	{
		cat("-----------------------------\n")
		info(sprintf("Checking %s", sig))
		genes = as.character(unlist(tab.ss2[, sig]))
		has.ensembl.and.isnt.in.10x = which((toupper(genes) %in% toupper(tab$Gene_Symbol.ss2)) & 
			(!toupper(genes) %in% toupper(tab$Gene_Symbol.10x)))
		
		can.replace = genes[has.ensembl.and.isnt.in.10x]
		replacement = as.character(tab$Gene_Symbol.10x[match(can.replace, tab$Gene_Symbol.ss2)])

		info(sprintf("Will replace:"))
		info(paste(can.replace, collapse=", "))
		info("Replacement:")
		info(paste(replacement, collapse=", "))
		cat("\n \n \n ")
		
		genes.clean = genes
		genes.clean[has.ensembl.and.isnt.in.10x] <- replacement
		length(genes.clean) = nrow(tab.ss2)
		tab.ss2[, sig]  = genes.clean
		#tab$Gene_Symbol.10x[match(psc.bad, tab$Gene_Symbol.ss2)]

	}
	return(tab.ss2)
}

make.ss2.10x.table <- function()
{
	ens2gene = read.delim("~/Desktop/projects/GutCircuits/csmillie/gene_sets/ensemble2gene.txt", header = F)
	colnames(ens2gene) = c("ENSEMBL", "Gene_Symbol.10x")
	ss2 = read.delim("~/ref/tables/rsem.txt")
	ss2.ensembl =  get.ensembl.for.ss2(as.character(ss2$gene_id))


	df = merge(ens2gene, ss2.ensembl, by="ENSEMBL") #, all=T)

	write.table(df, file="~/Desktop/10x_SS2_genes.txt", sep="\t", quote=F, row.names=F)
}

## Used to filter 10x signatures using smartseq 2 (SS2) plate data
## Use ENSEMBL IDs to map the genes to each other, and preserve rank using mean rank.
xref <- read.delim("~/ref/tables/biomart_xref.mm10.txt")
ucsc.rsem = read.delim("~/ref_temp/RSEM_mm10_ucsc_genomestudio_genes.genes.results")
ensembl.genes = read.delim("~/ref_temp/gene_id_table.txt") 	## mapping from 10x gene symbol to ENSEMBL id, 
															## extracted from /ahg/regevdata/users/ahaber/dev/cellranger/refdata-cellranger-1.0.0/mm10/genes/genes.gtf, 
															## using ~/dev/adam/rna_seq/rna_seq/gtf_to_id_table.py

get.consensus.sig <- function(sig.10x, sig.ss2, top.only.debug=0, all=F)
{

	if(top.only.debug>0)
	{
		sig.10x = sig.10x[1:top.only.debug]
		sig.ss2 = sig.ss2[1:top.only.debug]
	}

	#preprocess
	sig.10x = gsub(".", "-", sig.10x, fixed=T)
	sig.ss2 = gsub("Siglec5", "Siglecf", sig.ss2, fixed=T)

	# get ENSEMBL ids for 10x genes
	id.indices <- match(toupper(sig.10x), toupper(as.character(ensembl.genes$Gene_Symbol)))
	ENSEMBL = as.character(ensembl.genes$ENSEMBL_id)[id.indices]
	rank.10x = 1:length(sig.10x)
	df.10x = cbind.data.frame(sig.10x, ENSEMBL, rank.10x, stringsAsFactors=F)
	colnames(df.10x) = c("Gene_Symbol.10x", "ENSEMBL", "rank.10x")
	df.10x = df.10x[!is.na(df.10x$Gene_Symbol.10x),]
	

	extract.field = function(string,field=1,delim="_", fixed=T) {
		return(strsplit(string,delim, fixed=fixed)[[1]][field])
	}

	ambig.genes = grep("_", df.10x$Gene_Symbol.10x)
	if(length(ambig.genes)>0)
	{
		info("Parsing ambiguous gene names")
		genes_syms = unlist(lapply(as.character(df.10x$Gene_Symbol.10x)[ambig.genes], extract.field, field=1))
		ens_ids = unlist(lapply(as.character(df.10x$Gene_Symbol.10x)[ambig.genes], extract.field, field=2))
		df.10x$Gene_Symbol.10x[ambig.genes] <- genes_syms
		df.10x$ENSEMBL[ambig.genes] <- ens_ids
	}

	# get ENSEMBL ids for SS2 genes
	df.ss2 = get.ensembl.for.ss2(sig.ss2)

	

	# use ENSEMBL ids to match up genes
	
	# print("10X")
	# print(df.10x)

	# print("SS2")
	# print(df.ss2)

	df = merge(df.10x, df.ss2, by="ENSEMBL", all=all)
	#df = df[!is.na(df$ENSEMBL), ]
	
	# print(head(df))
	# print(head(df[, c("rank.ss2", "rank.10x")]))

	info("Sorting by combined rank")
	df$rank = rowMeans(as.matrix(df[, c("rank.ss2", "rank.10x")]))
	df = df[order(df$rank), ]

	info(sprintf("%s genes in common: ", nrow(df[!is.na(df$ENSEMBL),])))
	return(df)
}



get.ensembl.for.ss2 <- function(ss2.genes)
{
	df.ss2 = NULL
	rank = 1
	for (gene in ss2.genes)
	{
		#info(sprintf("Inspecting SS2 gene %s", gene))
		tab = xref[xref$Associated.Gene.Name == gene,]
		ens_id = as.character(unique(tab$Ensembl.Gene.ID))
		if(length(ens_id)==0)
		{
			ens_id = NA
			warn(sprintf("No ENSEMBL id for %s", gene))
		}else{
			if(length(ens_id) > 1){
				warn(sprintf("Multiple ENSEMBL ids for %s", gene))
				print(tab)
				#print(ens_id)
				### if the gene has multiple ENSEMBL ids, add them all to the table,
				### and then the merge will pick up only ENSEMBL ids that are found
				### in the 10x data, if any
				ens_ids = ens_id
				for( id in ens_ids)
				{
					if(is.null(df.ss2))
					{
						df.ss2 = data.frame(gene, id, rank)
						colnames(df.ss2) = c("Gene_Symbol.ss2", "ENSEMBL", "rank.ss2")
					}else{
						new.row = data.frame(t(c(gene, id, rank)))
						colnames(new.row) = c("Gene_Symbol.ss2", "ENSEMBL", "rank.ss2")
						df.ss2 = rbind.data.frame(df.ss2, new.row)
					}
				}
			}else{
				if(is.null(df.ss2))
				{
					df.ss2 = data.frame(gene, ens_id, rank)
					colnames(df.ss2) = c("Gene_Symbol.ss2", "ENSEMBL", "rank.ss2")
				}else{
					new.row = data.frame(t(c(gene, ens_id, rank)))
					colnames(new.row) = c("Gene_Symbol.ss2", "ENSEMBL", "rank.ss2")
					df.ss2 = rbind.data.frame(df.ss2, new.row)
				}
			}
		}
		rank = rank + 1
	}
	df.ss2$rank.ss2 = as.numeric(df.ss2$rank.ss2)
	df.ss2 = df.ss2[!is.na(df.ss2$Gene_Symbol.ss2),]
	return(df.ss2)
}



