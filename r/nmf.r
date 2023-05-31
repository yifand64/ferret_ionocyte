#
#	Run non negative matrix factorization to discover
#	basis components / metagenes that covary across samples.
#	http://nmf.r-forge.r-project.org/vignettes/NMF-vignette.pdf
#
run_nmf <- function(df, rank, output_name) {
  	require(NMF)
	res <- nmf(df, rank)
	layout(cbind(1, 2)
	basismap(res, subsetRow = TRUE)
	coefmap(res)
}
