# code from karthik.

library(Seurat)
table_file = "/Volumes/Boxcryptor/Dropbox/Postdoc/Projects/GutCircuits/Analysis/single_cell_miseq/tables/single_cell_counts_moshes_data_filtered.txt"
Counts_mat = read.table(table_file, header=TRUE, sep="\t", row.names=1, quote="")
rownames(Counts_mat) = as.character(Counts_mat[,1]); Counts_mat=Counts_mat[,c(-1,-2,-3)]
rownames(Counts_mat)=gsub("\"","", rownames(Counts_mat))

msh.gut =seurat(raw.data=log(Counts_mat+1),count.data=Counts_mat, ident.fxn=getStat3)
msh.gut=setup(msh.gut,project="MosheAdam",min.cells = 5,min.genes = 2000,calc.noise=FALSE,is.expr=0, min.counts=30)


msh.gut=mean.var.plot(msh.gut,y.cutoff = 1.5,y.high.cutoff = 40, x.low.cutoff = 0.5, fxn.x = mean, x.high.cutoff = 30, do.contour=FALSE)

msh.gut=pca(msh.gut,pcs.store = 20,pcs.print = 5,genes.print = 20)

msh.permutation = jackstraw.permutation.test(msh.gut, n.resamp=100, p.thres = 0.05, do.verbose=TRUE, seed.val=NULL)
msh.gut=run_tsne(msh.gut,pcs.use = 1:6,max_iter=1000,perplexity=20)
msh.gut=Mclust_dimension(msh.gut,1,2,reduction.use = "tsne",G.use = 2.5,set.ident = TRUE)


msh.gut.lgr = subsetData(msh.gut, cells.use=which.cells(msh.gut, c("Lgr5Hi", "Lgr5lo","Entero")))
msh.gut.lgr=mean.var.plot(msh.gut.lgr,y.cutoff = 1.6,y.high.cutoff = 40, x.low.cutoff = 0.5, fxn.x = mean, x.high.cutoff = 30, do.contour=FALSE)

msh.gut.lgr=pca(msh.gut.lgr,pcs.store = 20,pcs.print = 5,genes.print = 20)
msh.permutation = jackstraw.permutation.test(msh.gut.lgr, n.resamp=100, p.thres = 0.05, do.verbose=TRUE, seed.val=NULL)
