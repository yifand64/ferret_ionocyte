



# For some reason the metaprop function doesn't return the values in a nice way, this just runs the GLMM model and extracts
# estimate and confidence interval from the output
meta_prop <- function(x, n){
    library(metafor)
    library(meta)
    if(all(x==0)){return(list("x"=NA, "ci95.high"=NA, "ci95.low"=NA))}    
    if(any(n==0)){warn("Removing zero n samples (studies)"); non_zero_n = which(n>0); n = n[non_zero_n]; x = x[non_zero_n]}
    info(sprintf("Fitting.. [x=%s, n=%s]", paste(x, collapse=", "), paste(n, collapse=", ")))
    m = meta::metaprop(x, n, method = "GLMM")
    row.num = 5 #5 for random effects, 4 for fixed
    prop = as.numeric(strsplit(strsplit(capture.output(summary(m))[row.num], split=";", fixed=T)[[1]], split=" ")[[1]][8])
    ci95.low = as.numeric(strsplit(strsplit(capture.output(summary(m))[row.num], split=";", fixed=T)[[1]], split = "[", fixed=T)[[1]][2])
    ci95.high =  as.numeric(strsplit(strsplit(capture.output(summary(m))[row.num], split=";", fixed=T)[[1]], split = "]", fixed=T)[[2]][1])
    list("x"=prop, "ci95.high"=ci95.high, "ci95.low"=ci95.low)
}