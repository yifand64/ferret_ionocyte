# for a set of samples, calculate a score based on their expression of certain types of genes, 
# for example, cell-cycle genes.

library(plyr)
library(reshape2)



load.cluster.markers.into.gene.lists <- function(dir=".",top.n=50)
{
	marker.lists = list.files(dir)
	cat(sprintf("Reading %i marker lists..\n", length(marker.lists)))
	#cluster.names = unique(cluster.assignments)
	gene.lists = list()
	for(i in 1:length(marker.lists))
    {
        #markers.file = paste(cluster.names[i], paste("markers_",cluster.names[i], ".txt", sep=""), sep="/")
        markers.dir = marker.lists[i]
        index = grep("markers_", list.files(markers.dir))
        if(length(index) > 0)
        {
    		markers.list = list.files(markers.dir)[index]
    		markers.file = paste(markers.dir, markers.list, sep="/")
    		cat(sprintf("Reading markers for cluster %i of %i from %s ..\n", i, length(marker.lists), markers.file))
	        if (isTRUE(file.exists(markers.file)) & markers.list != "") 
	        {
	        	markers = read.delim(markers.file)
	        	genes = unlist(head(markers$GENE_SYMBOL, n = top.n))
        	}else{
        		cat("Shouldn't this never happen?\n")
        		cat(sprintf("WARN: No markers found in %s! \n" , markers.dir))
        		genes = c(NA)
        	}
        	# print("Loaded genes: ")
        	# print(genes)
        }else{
        	cat(sprintf("WARN: No markers found in %s! \n" , markers.dir))
        	genes = c(NA)
        }
        
        gene.lists[[i]] = genes
    }
    # make sure clusters with less than top.n markers dont cause an error
    print(gene.lists)
    lengths = lapply(gene.lists, function(x) length(x))
    print(lengths)
    n = lengths[[which.max(lengths)]]
    print("max length")
    print(n)

    # pad with NAs
    print("Padding with NAs")


	for(i in 1:length(marker.lists))
	{
		cat(sprintf("Markers list %i\n",i))
		
		if(length(gene.lists[[i]]) != n)
		{

			length(gene.lists[[i]]) = n
		}
		print(gene.lists[[i]])
		print("Length: ")
		print(length(gene.lists[[i]]))

		#print(cbind.data.frame(gene.lists[[i]], gene.lists[[i+1]]))
		#length(gene.lists[[i]]) = n
	}
	
	#convert to DF
	rval = as.data.frame(sapply(gene.lists, cbind.data.frame))

	#rval = data.frame(cbind(unlist(gene.lists)))
	
	colnames(rval) = marker.lists
	print(head(rval))
	return (rval)
}





