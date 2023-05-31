feature.plot.scale <- function(object,features.plot,pc.1=1,pc.2=2,cells.use=NULL,pt.size=1,cols.use=heat.colors(10),pch.use=16,reduction.use="tsne",nCol=NULL, show.labels=T) {
{
        
		if(is.null(cells.use))
		{
			cells.use = colnames(object@data)
		}
        dim.code="PC"
        if (is.null(nCol)) {
          nCol=2
          if (length(features.plot)>6) nCol=3
          if (length(features.plot)>9) nCol=4
        }         
        num.row=floor(length(features.plot)/nCol-1e-5)+1
        par(mfrow=c(num.row, nCol))
        if (reduction.use=="pca") {
          data.plot=object@pca.rot[cells.use,]
          dim.code="PC"
        }
        if (reduction.use=="tsne") {
          data.plot=object@tsne.rot[cells.use,]
          dim.code="tSNE_"
        }
        if (reduction.use=="ica") {
          data.plot=object@ica.rot[cells.use,]
          dim.code="IC"
        }
        
        ident.use=as.factor(object@ident[cells.use])
        data.plot$ident=ident.use
        x1=paste(dim.code,pc.1,sep=""); x2=paste(dim.code,pc.2,sep="")
        print(x1)
        print(data.plot)
        data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]
        data.plot$pt.size=pt.size
        data.use=data.frame(t(fetch.data(object,features.plot,cells.use = cells.use)))
        for(i in features.plot) {
			data.gene=na.omit(data.frame(data.use[i,]))
			data.cut=as.numeric(as.factor(cut(as.numeric(data.gene),breaks = length(cols.use))))
			data.col=rev(cols.use)[data.cut]
			par("mar"=c(5, 4, 4, 6))
			plot(data.plot$x,data.plot$y,col=data.col,cex=pt.size,pch=pch.use,main=i,xlab=x1,ylab=x2)
			data.range = c(min(data.gene), max(data.gene))
			if(max(data.gene > 5))
			{
				b=2
			}else{
				b=1
			}
			breaks = seq(0, floor(max(data.gene)), by=b)
			if(show.labels)
			{
				colorlegend(col=rev(cols.use), zlim=data.range, zval=breaks,  left=F, posy=c(0,0.8)) #main="Log2(count)",
			}else{
				colorlegend(col=rev(cols.use), zlim=NULL, zval=NULL, posy=c(0,0.8)) #main="Log2(count)",
			}
         }
        par(mfrow=c(1,1))
      }
}

# given a *NAMED* vector/dataframe of values (must be the same length as
# the number of cells) store it in the seurat object. useful for
# feature plots etc.
set.attribute <- function(seurat.obj, attrib.data)
{
    attrib.name=colnames(attrib.data)
    if(ncol(attrib.data) > 1)
    {
    	 seurat.obj@data.info = cbind(seurat.obj@data.info, attrib.data)
	}else{
		 seurat.obj@data.info[attrib.name]=attrib.data     
	}
   	
   	#print(colnames(seurat.obj@data.info))
    return (seurat.obj)
}

score.cells <- function(counts, gene.lists, alg="average", normalise=F, log =F)
{
	rval = NULL
	cat(sprintf("Scoring %i samples based on their expression of %i gene lists..\n", ncol(counts), length(gene.lists)))
	
	for (i in 1:length(gene.lists) ) {
		score.name = names(gene.lists)[i]
		cat(sprintf("Calculating %s score..\n", score.name))
		genes = gene.lists[[score.name]]
		print(genes)
		genes = gsub(" ", "", genes, fixed = TRUE) # god help you if there's spaces in the names of the genes
		counts.of.genes = counts[factor(rownames(counts)) %in% factor(genes), ]
		if(identical(alg,"average"))
		{
			scores = colMeans(counts.of.genes)
			if(normalise)
			{
				#min-max normalise [linear]
				#norm.scores = (scores - min(scores)) / (max(scores)-min(scores))
				
				if(log)
				{
					cat(sprintf("Taking log2..\n"))
					scores = log2(scores)
				}

				#sigmoid normalise
				norm.scores = scores / (1 + scores)
				scores = norm.scores
			}
		}else{
			stop(paste("Unsupported algorithm: " , alg))
		}
		if(is.null(rval))
		{
			rval = data.frame(scores)
			colnames(rval) = score.name
		}else{
			rval[score.name] = scores
		}
	}
	return (rval)
}