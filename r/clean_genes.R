clean_genes = function(genes){
  d  = data.frame(orig = genes, update=genes)
  annot = read.delim("~/Dropbox/HSPH/projects/ferret_ionocyte/hsph/analysis/0_gene_names/annot_table.txt")
  annot = annot %>% filter(query_gene %in% d$orig)
  d$update[match(annot$query_gene, d$update)] = annot$annot_geneonly
  d$update[is.na(d$update)] = d$orig[is.na(d$update)]
  d$update
}