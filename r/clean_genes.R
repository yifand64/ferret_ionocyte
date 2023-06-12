clean_genes = function(genes){
  d  = data.frame(orig = genes, update=genes)
  annot = read.delim("~/Desktop/haber/ferret/ferret_ionocyte/annot_table.txt")
  annot = annot %>% filter(query_gene %in% d$orig)
  d$update[match(annot$query_gene, d$update)] = annot$annot_geneonly
  d$update[is.na(d$update)] = d$orig[is.na(d$update)]
  d$update
}