

pheno.table.file = "~/ref/tables/HMD_HumanPhenotype_from_rebecca_jan_27.txt"

### convert the MSigDB environment from Livnat and Matan to mouse gene sets.
convert.go.env <- function(go.env, input.species='human')
{
    ids = names(go.env)
    n = length(ids)
    rval = list()
    i = 1
    for(id in ids)
    {
        gene.set = go.env[[id]]
        rval[[id]]  = gene.name.conversion(gene.set, input.species=input.species)
        # if(i>100)
        # {
        #   break
        # }
        i = i + 1
        info(sprintf("Gene-set %s of %s [%s%% complete]", i, n, 100*round(i/n, 3)))
    }
    rval <- clean.gos(rval) # remove GOs with too few or too many genes
    rval <- list2env(rval) # convert to an environment

    return(rval)
}


gene.name.conversion <- function(gene.set, input.species='mouse',return.homology.table=FALSE){
  
  # Carl's data
  library(mygene)
  library(reshape)

  Carl.gene.names <- read.delim('~/ref/tables/homologene_R68_human_mouse_121_orthologs_just_symbol.txt',header=T,stringsAsFactors = F)
  # Carl.gene.names = NULL
  # Carl.gene.names$Source <- 'Carl'
  
  # Mygene tool
  #info("Running mygene query")
  mygene.output <- queryMany(gene.set, scopes = 'symbol', species=c('mouse','human'))
  #print(mygene.output)	
  mygene.output <- subset(mygene.output, subset= (!is.na(taxid)))
  mygene.output <- subset(mygene.output, subset= (taxid == 9606 | (taxid == 10090 & query != symbol)))
  mygene.output <- unique(mygene.output[c('taxid','query','symbol')])
  #print(mygene.output)
  mygene.conversion <- na.omit(cast(mygene.output, query~taxid))
  #print(head(mygene.conversion))
  if(ncol(mygene.conversion) < 3)
  {
      warn("Lookup failed!")
      return(NA)
  }
  mygene.conversion <- mygene.conversion[, c(2,3)]
  colnames(mygene.conversion) <- c('symbolHS','symbolMM')
  # mygene.conversion$Source <- 'Mygene'
  
  # MGI homology table
  MGI.table <- read.delim(pheno.table.file, header=F,stringsAsFactors = F)
  MGI.table <- MGI.table[,c(1,4)]
  colnames(MGI.table) <- c('symbolHS','symbolMM')
  # MGI.table$Source <- 'MGI'
  
  #combining all three homology sources:
  if(!is.null(Carl.gene.names))
  {
		homology.gene.names <- unique(rbind(Carl.gene.names,mygene.conversion,MGI.table))
  }else
  {
  		homology.gene.names <- unique(rbind(mygene.conversion,MGI.table))
  }
  
  if (return.homology.table == TRUE){
    return(homology.gene.names)
  } else {
    if (input.species == 'mouse'){
      gene.set.converted <- homology.gene.names$symbolHS[homology.gene.names$symbolMM %in% gene.set]
    } else {
      gene.set.converted <- homology.gene.names$symbolMM[homology.gene.names$symbolHS %in% gene.set]
    }
    return(unname(gene.set.converted))
  }
}

load.scde.markers <- function(scde.markers.file)
{
	data <- read.delim(scde.markers.file,header=T)  
  	# path.split <- strsplit(scde.markers.file, '\\/|_markers.txt')[[1]]
  	exp.name = strsplit(scde.markers.file, "markers_")[[1]][2]
  	exp.name = gsub(".txt", "", exp.name)
  	info(sprintf("Loaded markers for %s", exp.name))
  	return (data)
}



GSEA.prepranked.from.ranked.table <- function(data, exp.name, homology.table=NULL){
  
  if(is.null(homology.table))
  {
  		info("No homology table provided, converting gene names to human for GSEA")
  		curr.homology.table <- gene.name.conversion(data$GENE_SYMBOL, input.species='mouse',return.homology.table=TRUE)
  }else
  {
  		info("Using provided homology table")
  		curr.homology.table = homology.table
  }
	
  rescale = Vectorize(function(x){x/(x+1)})
  rescale.minmax = function(x)(x-min(x))/(max(x)-min(x))

  info("Matching gene symbols to homology table")
  data$Rank = 1:nrow(data)
  
  #  map ranks to 0-10 so highest ranked gene has score 10, lowest 0
  data$Rank.Score = 10*rescale.minmax(1/data$Rank)

  data <- merge(x = data,y = curr.homology.table,all.x=T,by.x='GENE_SYMBOL',by.y='symbolMM')
  data = data[order(data$Rank), ]

  #print(head(na.omit(data[c('GENE_SYMBOL', 'symbolHS','Rank')]), n=20))
  output.file = paste0(exp.name,'_Preranked.rnk')
  info(sprintf("Writing ranked list to %s", output.file))

  

  ranked.list <- na.omit(data[c('symbolHS','Rank.Score')])
  
  write.table(ranked.list, 
  	file = output.file,
  	col.names=F,row.names=F,quote=F,sep='\t')
}

