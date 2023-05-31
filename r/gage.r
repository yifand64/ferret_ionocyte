	
library(gage)
library(gageData)
data(kegg.gs) #initalizes Kegg library; other possible selections:
# data(kegg.sets.hs)
# data(go.sets.hs)
# data(carta.hs)
# data(kegg.sets.mm)
# data(go.sets.mm)
# data(kegg.sets.rn)
# data(go.sets.rn)
# data(kegg.sets.sc)
# data(go.sets.sc)
# data(sigmet.idx.hs)
# data(go.subs.hs)
# data(sigmet.idx.mm)
# data(go.subs.mm)
# data(sigmet.idx.rn)
# data(go.subs.rn)
# data(sigmet.idx.sc)
# data(go.subs.sc)


data(egSymb) #initalizes Entrez ID -> gene symbol mappings 

# run_gage <- function(df, ref_cluster_id) { #non ref id automatically assinged to sample cluster
  
# }



# ref.columns is the indices of the 'sample group'
# all others are considered background.
run_gage <- function(counts.file, ref.columns, min.q.val=0.1)
{
	df <- readExpData(counts.file, row.name =1)
	
	rownames(df) <- toupper(rownames(df)) #capitalize all gene IDs for proper mapping for KEGG
	kegg.gs.sym <- lapply(kegg.gs, eg2sym)
	kegg.p <- gage(df, gsets = kegg.gs.sym, ref = ref.columns, compare = 'unpaired')
	kegg.p.sig_up <- kegg.p$greater[, "q.val"] < min.q.val & !is.na(kegg.p$greater[, "q.val"])
	kegg.p.sig_down <- kegg.p$lesser[, "q.val"] < min.q.val  & !is.na(kegg.p$lesser[, "q.val"])
	
	cat(sprintf("%s pathways upregulated [q < %s]", nrow(kegg.p), min.q.val))
	# head(kegg_p$greater) #pathway groups up in sample group
	# head(kegg_p$less) #pathway groups down in sample group
	return (kegg.p)
}