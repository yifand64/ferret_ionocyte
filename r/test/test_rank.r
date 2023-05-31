source("~/dev/adam/rna_seq/r/util.r")
script("marker_compare")
library(readr)

cell.type = "Quiescent.Stem"
background = c("Cycling.Stem", "Goblet")

tissue = "SI"
cell.markers.file = paste("~/dev/adam/rna_seq/r/gui/online/Gut\ Atlas\ GUI/resources/", 
							tissue,"_markers/markers_", cell.type, "_full.txt", sep="")
markers <- read_tsv(cell.markers.file)

load("~/dev/adam/rna_seq/r/gui/online/Gut\ Atlas\ GUI/resources/SI_780.Rdata")
norm.counts = seurat.obj@data
clusters = seurat.obj@data.info$DBclust.ident


info(sprintf("Comparing against %s", background))  
keep.cols = c("GENE_SYMBOL", 
    paste("Conservative.Estimate", cell.type, "vs", background, sep="_"), 
    paste("p.adj", cell.type, "vs", background, sep="_"),
    paste("Corrected.Z.score", cell.type, "vs", background, sep="_"))
print(keep.cols)
markers = data.frame(markers[ , keep.cols])
print(dim(markers))

if(length(background) > 1)
{
    info("Ranking..")
    markers = do.rank(markers, cell.type, clusters, norm.counts)
    print(head(markers))
}else
{
	print(head(markers))
}

print(concise.cols[!concise.cols %in% colnames(markers)])
write.table(markers[, concise.cols], file="markers_test.txt", sep="\t", quote=F, row.names=F)