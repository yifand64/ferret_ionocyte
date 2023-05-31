library(Hmisc)
library(ggplot2)
library(mgcv)
library(viridis)
library(parallel)
library(foreach)
library(doParallel)
library("ElPiGraph.R")
library(igraph)
library(future.apply)
library(progressr)


### slightly adjust the previous the monocle-based function to work with any pseudotime ordering.
### Draw a smoothed (loess) heatmap of the given genes across the given cells.
## pt - pseudotime
transition.heatmap <- function(obj, cells.use, pt, transition.name, genes,
	reverse=F,
	pdf.width=7, 
	pdf.height=5.5, 
	colors=rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)), 
	disp.range=3, 
	num.expr.bins=30, ## order genes by pseudo-timepoint where expression is highest
	loess.span=0.5, 
	cex.row = 3,
	pseudotime.cols=cubehelix1.16)
{
	d = GetAssayData(obj)[genes, cells.use]
	keep_genes = 1:nrow(d)
	if(any(rowSums(d)==0)){
		n_zero=sum(rowSums(d)==0); 
		keep_genes = which(rowSums(d) > 0)
		d = d[keep_genes,]; 
		warn(sprintf("Removed %s zero rows", n_zero))
	}
	reorder = if(reverse) rev(order(pt)) else order(pt)
	du = d[, reorder]
	ac = data.frame(annCol = pt[reorder]); colnames(ac) = "pseudo-time" #set up coloring to indicate pseudotime
	if(nrow(du) < 2){stop("Less than 2 significant hits. Stopping")}	
	
	genes_per_bin = round(length(pt) / num.expr.bins)
	bins = list(cut2(pt[reorder], m=genes_per_bin)) #list(cut(pt[reorder], breaks = num.expr.bins))
	du = smooth.data(du, loess.span=loess.span) 
	du = t(scale(t(du)))
	## order genes by pseudo-timepoint where expression is highest (split pseudotime interval into bins. which bin has the highest mean expression?)
	g = data.frame(gene=rownames(du), gene_order=unlist(apply(du, 1, function(x) which.max(aggregate(unlist(x), by = bins, FUN=median)$x))))
	rownames(g) = rownames(du)
	du = if(reverse) du[rev(order(g$gene_order)),] else du[order(g$gene_order),]
	g = g[rownames(du),]
	g$gene = NULL
	## clip values 
	du[du > disp.range] = disp.range
	du[du< -disp.range] = -disp.range

	pseudo.cols = list("pseudo-time"=if(reverse) rev(pseudotime.cols) else pseudotime.cols)

	info(sprintf("Drawing transition [%s] heatmap for %s genes, %s cells", transition.name ,nrow(du), ncol(du)))
	write.table(rownames(du), file=sprintf("%s_genes_ordered.txt", transition.name), sep="\t", quote=F)
	aheatmap(du, filename=sprintf("%s_autocor.pdf", transition.name), 
		main=sprintf("%s [%s cells]", transition.name, length(cells.use)),
		width=pdf.width, height=pdf.height, Colv = NA, Rowv = NA, annCol = ac, annColors = pseudo.cols, 
		color = colors, cexRow = cex.row, annRow = g, cexCol = 0)
}


pt_plot = function(obj, gene, cell, cell.ids=NULL, vln=F, vln_bins=10, pt_name="pseudotime_paga")
{
	info(sprintf("Pseudotime plot. Cell=%s, Gene=%s", cell, gene))
	blank_x_axis = theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  if(is.null(cell.ids)){
    #cell.ids = obj$PGT_branch_ID_named %in% c("Precursor", cell)
    cell.ids = 1:ncol(GetAssayData(obj))
  }
  cells.use = colnames(obj)[cell.ids]
  pseudotime = obj@meta.data[,pt_name][cell.ids]
  expr = GetAssayData(obj)[gene, cell.ids]
  x = data.frame(expr, pseudotime)
  cb = stat_summary(geom = "crossbar", show.legend=FALSE, width=0.75, fatten=1, color="black",
					fun.data = function(x){ return(c(y=mean(x), ymin=mean(x), ymax=mean(x))) })
  if(vln){
    			x$pseudotime_bin = cut2(x$pseudotime, g = vln_bins)
    			ggplot(x, aes(x=pseudotime_bin, y=expr, fill=pseudotime_bin)) + 
    					geom_violin(scale="width",adjust=1,trim=TRUE, alpha=0.6, color=NA) + stat_summary(color="black", size=0.25) + scale_fill_manual("", values=colorRampPalette(matplotlib.viridis)(vln_bins), guide=F) +
    					geom_quasirandom(varwidth = TRUE, size=0.5, stroke=0.1, colour="grey50", show.legend=FALSE, method = "pseudorandom")  + 
    					stat_summary(color="black", fun.data="mean_cl_boot", fun.args=list(conf.int=.99), size=0.8, geom="errorbar", width=0.1) + #cb + 
    					geom_smooth(method="gam", formula=y ~ s(x, bs = "cs", k=4)) + ylab(gene) + blank_x_axis + xlab(sprintf("Binned Pseudotime (%s bins)", vln_bins))
    }else{
    			ggplot(x, aes(x=pseudotime, y=expr)) + geom_point(stroke=0.1, size=0.5) + 
    						ylab(gene) + 
    						geom_smooth(method="gam", formula=y ~ s(x, bs = "cs", k=4), size=2) 
  }
}


pt_multi_plot = function(obj, genes, cell, pt_name="pseudotime_paga")
{
	info(sprintf("Pseudotime plot. Cell=%s, Genes=%s", cell, genes))
    cell.ids = obj$PGT_branch_ID_named %in% c("Precursor", cell)
    cells.use = colnames(GetAssayData(obj))[cell.ids]
    x = FetchData(obj, c(genes, pt_name))[cell.ids, ]
  	x = melt(x, id.vars=pt_name)
	ggplot(x, aes_string(x=pt_name, y="value", color="variable")) + geom_point(stroke=0.1, size=0.5) + 
        scale_color_manual("", values=default.cols(length(genes))) + 
				ylab("Expression level, log2(TPM+1)") + xlab(sprintf("Pseudotime (Basal to %s)", cell)) + 
        geom_smooth(method="gam", formula=y ~ s(x, bs = "cs", k=4), size=1) 
}


# One method for detecting when a gene turns on is when LM estimate (or lower bound) goes above zero.
# watch out for giant file sizes.
# http://www.win-vector.com/blog/2014/05/trimming-the-fat-from-glm-models-in-r/
#  As many R users know (but often forget), a glm model object carries a copy of its training data by default. You can use the settings y=FALSE and model=FALSE to turn this off.
fit_gam = function(x = NULL, y = NULL, nfit_pts = 100, plot=T, print.dot=F, verbose=T, gene_name="Gene", knots=5)
{
   tryCatch({   
              if(verbose){info(sprintf("[fit_gam] Fitting %s", gene_name))}
              if(is.null(x)){x= seq(0, pi * 2, 0.1)}
              if(is.null(y)){y=sin(x) + rnorm(n = length(x), mean = 0, sd = 1)}
              gam_y = mgcv::gam(as.formula(sprintf("y ~ s(x, k=%s)", knots)), method = "REML")
              s = summary(gam_y)
              p_value = summary(gam_y)$s.pv # approximate p-value for smooth terms.
              x_new = seq(min(x), max(x), length.out = nfit_pts) 
              p = predict(gam_y, data.frame(x=x_new), se.fit = TRUE)
              upr = p$fit + (3 * p$se.fit) # 99% CI
              lwr = p$fit - (3 * p$se.fit) # 99% CI
              m = data.frame(x_new, p$fit, upr, lwr)
               if(plot) {g = ggplot(m, aes(x=x_new, y=p.fit)) + geom_ribbon(aes(ymin = lwr, ymax = upr ), fill = "grey70") + 
                            ggtitle(sprintf("GAM. p: %s, Dev. expl: %s, R-sqr: %s", signif(s$s.pv, 3), signif(s$dev.expl, 3), signif(s$r.sq, 3))) + 
                            ylab(sprintf("%s expression, log2(TPM+1)", gene_name)) + xlab("Pseudotime") + cowplot::theme_cowplot()  +
                            geom_point(data=data.frame(x,y), aes(x=x, y=y)) + geom_line(color="red", size=0.75) } else{ g = NA}
              if(print.dot) {cat(".")}
              list("plot"=g, "fit"=m, "p.value"=p_value, "gene"=gene_name, "r_sqr"=s$r.sq, "dev_expl"=s$dev.expl)}
        , interrupt = function(ex) {
          if(verbose){
                warn("[fit_gam] An interrupt was detected.")
                warn(ex)
          }
          list("plot"=NA, "fit"=NA, "p.value"=NA, "gene"=gene_name, "r_sqr"=NA, "dev_expl"=NA)
        }, error = function(ex) {
          if(verbose){
                error("[fit_gam] An error was detected.")
                error(ex)
          }
          list("plot"=NA, "fit"=NA, "p.value"=NA, "gene"=gene_name, "r_sqr"=NA, "dev_expl"=NA)
        }, finally = {})
}


fit_all = function(obj, 
                  genes, 
                  model = "gam",
                  pseudotime_id, 
                  cell.ids=NULL, 
                  nfit_pts = 100, 
                  n.cores=1, 
                  verbose=F, 
                  do_plots=F, 
                  knots=-1, 
                  normalize_pt=T)
{
      start_time <- Sys.time()
      if(!model %in% c("gam", "lm")){stop("Model must be one of lm or gam")}
      if(do_plots){warn("Having do_plots on will dramatically slow things down for high numbers of genes")}
      info("Fetching data from seurat object")
      if(is.null(cell.ids)){cell.ids=colnames(obj)}
      x = FetchData(obj, c(genes, pseudotime_id))[cell.ids, ] #"PAGA_1", "PAGA_2"
      pt = x[, pseudotime_id]
      if(normalize_pt){
        info("Normalizing pseudotime to the interval 0-1")
        pt = range01(pt)
      }
      # return values:
      plots = list()
      fits = list()
      models = list()
      nelements = length(genes)
      r_sqr    = rep(NA, nelements)
      pval     = rep(NA, nelements)
      dev_expl  = rep(NA, nelements)
      fit = matrix(nrow=nfit_pts, ncol=length(genes)) # transpose things bc iterative row access is slow.
      colnames(fit) = genes
      info(sprintf("Fitting %s %s models using %s cores", length(genes), model, n.cores))
      if(n.cores <= 1){
              for(i in 1:length(genes)) 
              {
                   g = genes[i]
                   done = scales::percent(i/length(genes), accuracy=0.1)
                   info(sprintf("Fitting gene %s,  %s of %s (%s)", g, i, length(genes), done))
                   expression = x[, g]
                   if(model == "gam"){
                      y = fit_gam(x=pt, y=expression, nfit_pts = nfit_pts, plot=do_plots, verbose=verbose, gene_name=g, knots=knots)
                   }else{
                      y = fit_lm(x=pt, y=expression, nfit_pts = nfit_pts, plot=do_plots, verbose=verbose, gene_name=g)
                   }
                   if(do_plots) {plots[[g]] = y$plot}
                   fits[[g]] = y$fit
                   genes[i] = y$gene
                   r_sqr[i] = y$r_sqr
                   pval[i]  = y$p.value
                   dev_expl[i]=y$dev_expl
                   if(!is.null(y$fit)){fit[,g] =  y$fit$lwr}else{fit[, g] = NA}
              }
              rv = list("fit"=fit, "plots"=plots, "genes"=genes, "r_sqr"=r_sqr, "p.value"=pval)
      }else{
              info(sprintf("Setting up parallelization stuff (n.cores=%s)", n.cores))
              export.obj.names = c("info", "warn", "error", "flog.info", "flog.warn", "flog.error", "fit_gam", "fit_lm")
              cl= makePSOCKcluster(n.cores, outfile="", setup_strategy = "sequential")
              registerDoParallel(cl, n.cores)
              d = t(x[, genes])
              rval=tryCatch({   
                  chunksize = nrow(d)/n.cores ### break the list of genes into chunks, one for each thread.
                  vals = split(1:nrow(d), ceiling(seq_along(1:nrow(d))/chunksize))
                  p <- foreach(run.id=1:n.cores,.export=export.obj.names, .combine=c) %dopar% {
                    v = vals[[run.id]]
                    init.val = v[1]
                    lapply(v, function(i) {
                          library(ggplot2)
                          gene.num = i - init.val
                          expression = unlist(d[i,])
                          if(round(length(v)/10) > 1){
                          if(gene.num %% round(length(v)/10) ==0  && run.id==1) # only the first chunk (i.e first process) posts updates
                          {
                              cat("\n")
                              info(sprintf("%s%% complete. Tested %s [%s of %s for this core]", 
                              100*round(gene.num/length(v), 3), rownames(d)[i], i, length(v)))
                          }}
                          if(model == "gam"){
                              y =fit_gam(x=pt, y=expression, nfit_pts = nfit_pts, plot=do_plots, verbose=verbose, print.dot=T, gene_name=rownames(d)[i], knots=knots)
                          }else{
                              y= fit_lm(x=pt, y=expression, nfit_pts = nfit_pts, plot=do_plots, verbose=verbose, print.dot=T, gene_name=rownames(d)[i])
                          }
                          y
                        })
                  }
                  list("results"=p, "cluster"=cl)
              }
              , interrupt = function(ex) {
                cat("An interrupt was detected.\n");
                print(ex);
                list("cluster"=cl)
              }, error = function(ex) {
                cat("An error was detected.\n");
                print(ex);
                list("cluster"=cl)
              }, finally = {
                list("cluster"=cl)
              })           

              cat("\n") #important for a new line after all the dots above
              nelements = length(rval$results)
              info(sprintf("Processing results for %s genes", nelements))
              genes     = rep(NA, nelements)
              plots     = as.list(rep(NA, nelements))
              fit   = as.list(rep(NA, nelements))
              for(i in 1:nelements) {
                  result = rval$results[[i]]
                  if(length(result$fit) > 1){fit[[i]] = result$fit$p.fit}else{
                      fit[[i]] = rep(NA, nfit_pts)
                  } # if length is 1, it's NA because the fit failed
                  if(do_plots) {plots[[i]] = result$plot}
                  genes[i] = result$gene
                  r_sqr[i] = result$r_sqr
                  pval[i]  = result$p.value
                  dev_expl[i]=result$dev_expl
              }
              fit   = data.frame(do.call(cbind, fit))
              colnames(fit) = genes
              info("Unregistering parallel backend and closing socket connections")
              stopCluster(rval$cluster)
              registerDoSEQ()
              closeAllConnections()
              info("Done")
              rv = list("fit"=fit, "plots"=plots, "genes"=genes, "r_sqr"=r_sqr, "dev_expl"=dev_expl, "p.value"=pval)
      }
      end_time <- Sys.time()
      info(sprintf("Computed GAM fits for %s genes in %s", length(genes), lubridate::as.duration(end_time - start_time)))
      return(rv)
}



# One method for detecting when a gene turns on is when lower bound goes above zero.
# watch out for giant file sizes.
# http://www.win-vector.com/blog/2014/05/trimming-the-fat-from-glm-models-in-r/
#  As many R users know (but often forget), a glm model object carries a copy of its training data by default. You can use the settings y=FALSE and model=FALSE to turn this off.
fit_lm = function(x = NULL, 
                  y = NULL, 
                  nfit_pts = 100, 
                  plot=T, 
                  print.dot=F, 
                  verbose=T, 
                  quantile_reg=F,
                  gene_name="Gene")
{
   tryCatch({   
              if(is.null(x)){x = seq(0, pi * 2, 0.1)}
              if(is.null(y)){y = 1-0.05*x^2.5 + rnorm(n = length(x), mean = 0, sd = 1)}
              x_new  = seq(min(x), max(x), length.out = nfit_pts) 
              if(quantile_reg){ 
                  model = quantreg::rq(y ~ poly(x, 2, raw=TRUE), tau = .5, ci=T)
                  s = quantreg::summary.rq(model, se = "boot")$coefficient
                  pv_lin    = s["poly(x, 2, raw = TRUE)1", "Pr(>|t|)"]
                  pv_sqr    = s["poly(x, 2, raw = TRUE)2", "Pr(>|t|)"]
                  fit = predict(model, data.frame(x=x_new), se.fit = TRUE)
                  r_sqr = NA
                  estimate = fit[,1]
                  upr      = fit[,2]
                  lwr      = fit[,3]
              }else{
                  model = lm(y ~ poly(x, 2, raw=TRUE))
                  s = summary(model)
                  p = s$coefficients[,4]
                  pv_lin = p[2]
                  pv_sqr = p[3]
                  fit = predict(model, data.frame(x=x_new), se.fit = TRUE)
                  estimate = fit$fit
                  upr      = fit$fit + (2 * fit$se.fit) # 95% CI
                  lwr      = fit$fit - (2 * fit$se.fit) # 95% CI
                  r_sqr = s$r.squared
              }
              title    = sprintf("LM x^2 fit \n Linear p: %s, Quadratic p: %s, Multiple R-sqr: %s", 
                signif(pv_lin, 3), signif(pv_sqr, 3), r_sqr)
              m = data.frame(x_new, estimate, upr, lwr)
               if(plot) {g = ggplot(m, aes(x=x_new, y=estimate)) + geom_ribbon(aes(ymin = lwr, ymax = upr ), fill = "grey70") + 
                            ggtitle(title) + ylab(sprintf("%s expression, log2(TPM+1)", gene_name)) + xlab("Pseudotime") + cowplot::theme_cowplot()  +
                            geom_point(data=data.frame(x,y), aes(x=x, y=y)) + geom_line(color="red", size=0.75)} else{ g = NULL}
              if(print.dot) {cat(".")}
              list("plot"=g, "fit"=m, "p_val_linear"=pv_lin, "p_val_quadratic"=pv_sqr, "gene"=gene_name, "r_sqr"=r_sqr, "dev_expl"=NA)}
        , interrupt = function(ex) {
          if(verbose){
                futile.logger::flog.warn("[fit_lm] An interrupt was detected.")
                futile.logger::flog.warn(ex)
          }
          list("plot"=NULL, "fit"=NULL, "p_val_linear"=NULL, "p_val_quadratic"=NULL, "gene"=gene_name, "r_sqr"=NA, "dev_expl"=NA)
        }, error = function(ex) {
          if(verbose){
                futile.logger::flog.error("[fit_lm] An error was detected.")
                futile.logger::flog.error(ex)
          }
          list("plot"=NULL, "fit"=NULL, "p_val_linear"=NULL, "p_val_quadratic"=NULL, "gene"=gene_name, "r_sqr"=NA, "dev_expl"=NA)
        }, finally = {})
}

# Folded this functionality into fit_all_gams
# fit_all_lms = function(obj, genes, pseudotime_id, cell.ids=NULL, nfit_pts = 100, n.cores=1, verbose=F, do_plots=F, normalize_pt=T, use_lower_bound=T, quantile_reg=F)
# {
#       start_time <- Sys.time()
#       info("Fetching data from seurat object")
#       if(is.null(cell.ids)){cell.ids=colnames(obj)}
#       x = FetchData(obj, c(genes, pseudotime_id))[cell.ids, ] #"PAGA_1", "PAGA_2"
#       pt = x[, pseudotime_id]
#       if(normalize_pt){
#         info("Normalizing pseudotime to the interval 0-1")
#         pt = range01(pt)
#       }
#       if(use_lower_bound){info("use_lower_bound is ON. Will return the fit for the lower end of the 95-CI")}
#       if(quantile_reg){info("quantile_reg is ON. Will return the quantile regression fit for the median")}
#       # return values:
#       plots = list()
#       fits = list()
#       models = list()
#       lm_fit = matrix(nrow=nfit_pts, ncol=length(genes)) # transpose things bc iterative row access is slow.
#       colnames(lm_fit) = genes
#       info(sprintf("Fitting %s LMs using %s cores", length(genes), n.cores))
#       if(n.cores <= 1){
#               for(i in 1:length(genes)) 
#               {
#                  g = genes[i]
#                  done = scales::percent(i/length(genes), accuracy=0.1)
#                  info(sprintf("Fitting gene %s,  %s of %s (%s)", g, i, length(genes), done))
#                  expression = x[, g]
#                  y = fit_lm(x=pt, y=expression, nfit_pts = nfit_pts, plot=do_plots, verbose=verbose, gene_name=g, quantile_reg=quantile_reg)
#                  if(do_plots) {plots[[g]] = y$plot}
#                  models[[g]] = y$model
#                  fits[[g]] = y$fit
#                  if(use_lower_bound){
#                     if(!is.null(y$fit)){lm_fit[,g] =  y$fit$lwr}else{lm_fit[, g] = NA} # could also use the lower bound
#                  }else{
#                     if(!is.null(y$fit)){lm_fit[,g] =  y$fit$estimate}else{lm_fit[, g] = NA}
#                  }
#               }
#               end_time <- Sys.time()
#               info("Extracting R-squared values")
#               r_sqr = lapply(models, function(x) summary(x)$adj.r.squared)
#               info(sprintf("Computed LM fits for %s genes in %s seconds", length(genes), signif(end_time - start_time, 2)))
#               return(list("lm_fit"=lm_fit, "plots"=plots, "models"=models, "genes"=genes, "r_sqr"=r_sqr))
#       }else{
#               warn("UNTESTED - unsure if it works..")
#               info(sprintf("Setting up parallelization stuff (n.cores=%s)", n.cores))
#               export.obj.names = c("info", "warn", "error", "flog.info", "flog.warn", "flog.error", "fit_lm" )
#               cl<-makePSOCKcluster(n.cores, outfile="", setup_strategy = "sequential")
#               d = t(x[, genes])
#               registerDoParallel(cl, n.cores)
#               rval=tryCatch({   
#                   chunksize = nrow(d)/n.cores ### break the list of genes into chunks, one for each thread.
#                   vals = split(1:nrow(d), ceiling(seq_along(1:nrow(d))/chunksize))
#                   p <- foreach(run.id=1:n.cores,.export=export.obj.names, .combine=c) %dopar% {
#                     v = vals[[run.id]]
#                     init.val = v[1]
#                     lapply(v, function(i) {
#                           library(ggplot2)
#                           gene.num = i - init.val
#                           expression = unlist(d[i,])
#                           if(gene.num %% 200==0  && run.id==1)
#                           {
#                             info(sprintf("%s%% complete. Tested %s [%s of %s]", 
#                               100*round(gene.num/length(v), 3), rownames(d)[i], i, length(v)))
#                           }
#                           fit_lm(x=pt, y=expression, nfit_pts = nfit_pts, plot=do_plots, print.dot=T, verbose=verbose, gene_name=rownames(d)[i], quantile_reg=quantile_reg)
#                         })
#                   }
                  
#                   #rv = data.frame(matrix(unlist(p), nrow=nrow(d), byrow=T),stringsAsFactors=FALSE)
#                   #colnames(rv) = names(p[[1]])  
#                   list("results"=p, "cluster"=cl)
#               }
#               , interrupt = function(ex) {
#                 cat("An interrupt was detected.\n");
#                 print(ex);
#                 list("cluster"=cl)
#               }, error = function(ex) {
#                 cat("An error was detected.\n");
#                 print(ex);
#                 list("cluster"=cl)
#               }, finally = {
#                 list("cluster"=cl)
#               })              
#               cat("\n")
#               info("Unregistering parallel backend")
#               stopCluster(rval$cluster)
#               registerDoSEQ()
#               info("Done")
#       }
#       #lm_fit[lm_fit < 0] = 0
#       #lm_fit[lm_fit > max] = max

#       end_time <- Sys.time()
#       info(sprintf("Computed LM fits for %s genes in %s seconds", length(genes), signif(end_time - start_time, 2)))
#       return(rval$results)
# }

# wrote this for the neurogenesis paper with Lora.
ordered_heatmap = function(m, nonzero_thresh = 0.25, 
                              rsqr_thresh = 0.025, 
                              filter_to_gene_list=NULL, 
                              scale_expression=T, 
                              cols= NULL,
                              width = 10,
                              height = 10,
                              remove_initially_on=T, # remove genes that are on at the start
                              show_rownames = T,
                              verbose=T,
                              fontsize_row=4,
                              clip = 3, # maximum value to improve color contrast
                              rank_by = "first_nonzero",
                              destfile="pseudotime_NG_TFs_order.pdf")
{
    d = data.frame(m$genes, unlist(m$r_sqr)) 
    rownames(d) = m$genes
    colnames(d) = c("gene", "r_sqr")
    info(sprintf("Got fit info for %s genes", nrow(d)))
    if(!all(rank_by %in% c("max_loc", "max_size", "first_nonzero"))){stop("rank_by must be one or more of: max_loc, max_size, first_nonzero")}
    if(!is.null(rsqr_thresh)) {d = subset(d, r_sqr > rsqr_thresh); info(sprintf("   %s genes pass R-squared cut-off [%s]", nrow(d), rsqr_thresh))}
    if(!is.null(filter_to_gene_list)) {d = d[d$gene %in% filter_to_gene_list,]; info(sprintf("   %s genes are in provided filter list [%s genes]", nrow(d), length(filter_to_gene_list)))}
    
    df = m$gam_fit[, match(d$gene, m$genes)] # get the fitted estimates from our LMs (or GAMs, but not implemented atm)
    
    if(!is.null(nonzero_thresh)) {
      nonzero_filter = apply(df, 2, function(x) any(x > nonzero_thresh))
      d = d[nonzero_filter,]
      df = df[, nonzero_filter]
      info(sprintf("   %s genes pass nonzero expression cut-off [%s]", nrow(d), nonzero_thresh))
    }

    d$gene_order = NA
    for(i in 1:nrow(d)){
      #print(d$gene[i])  
      d$gene_order[i] = min(which(df[, i] > nonzero_thresh)) # this is the first point in the x_new vector where the lower bound is above a threshold
    }

    df = t(df)
    ord = data.frame(apply(df, 1, which.max),  # order first by location of max, and then by size of max, and by location of first non-zero value
                     apply(df, 1, max), 
                     apply(df, 1, function(x) min(which(x > nonzero_thresh)))) 
    colnames(ord) = c("max_loc", "max_size", "first_nonzero")
    rownames(ord) = rownames(df)

    if(remove_initially_on){
      ord = subset(ord, first_nonzero > 1); 
      df = df[rownames(ord),]; 
      info(sprintf("   %s genes are not initially on", nrow(df)))
    }
    
    
    if(verbose){ggsave(ggplot(ord, aes(x=max_loc, y=first_nonzero)) + geom_point() + stat_density_2d() , filename="order_metrics.pdf")} # DEBUG plot

    if(scale_expression){ df = t(scale(t(df)))}
    
    if(clip > 0){
      df[df > clip] = clip
      df[df < -clip] = -clip
    }

    #iloc = ord[, "max_loc"]
    # flip_sort_rows = iloc < (nrow(df)/2) # this part seemed to be important for LM but not for GAMs? have to check.
    # ord$max_size[flip_sort_rows] = -1 * ord$max_size[flip_sort_rows] # flip the sort for one side to get the ordering right
    
    #ord = ord[order(ord$first_nonzero, ord$max_size),]
    ord = ord[order(ord[,rank_by]),]
    annot_row = data.frame(d[rownames(ord),]$r_sqr)
    colnames(annot_row) = "R-sqr"
    rownames(annot_row) = rownames(ord)
    if(is.null(cols)){if(scale_expression) {cols=colorRampPalette(rev(brewer.pal(11, "RdBu")))(25)}else{
      cols = image_cols
    }}
    info(sprintf("   saving heatmap to %s", destfile))
    #if(nrow(ord) > 500){show_rownames=F}
    df = df[rownames(ord),]
    pheatmap(df, cluster_cols = F, cluster_rows = F, 
      annotation_row = annot_row, color = cols, show_rownames=show_rownames,
      filename=destfile, width=width, height=height, fontsize_row = fontsize_row, border_color = NA)
    write.table(data.frame("rank"=1:length(rownames(ord)), "gene"=rownames(ord)), file=gsub(".pdf", ".txt", destfile), sep="\t", quote=F, row.names=F)
}


### use autocorrelation to find genes significantly associated with each edge.
get.edge.genes <- function(obj, cells.use, pt, transition.name, 
  reverse=F,
  lag=25, 
  min.expr=1, 
  min.cells=5, 
  n.genes=100, 
  n.genes.tfs=20, 
  B=100, 
  show.p=0.05, 
  n.cores=1)
{
  d = GetAssayData(obj)[, cells.use]
  reorder = if(reverse) rev(order(pt)) else order(pt)
  d = d[, reorder]
  ac = data.frame(annCol = pt[reorder]); colnames(ac) = "pseudo-time" #set up coloring to indicate pseudotime
  d.use = d[Matrix::rowSums(d>min.expr) > min.cells,]
  info(sprintf("Calculating autocorrelations and p-values for %s genes, %s cells [B=%s]", nrow(d.use), ncol(d.use), B))

  cors = calc.autocors(d.use, B=B, n.cores=n.cores)
  cors = na.omit(cors)
  cors = cors[order(cors$ac, decreasing=T),]
  cors$fdr = p.adjust(cors$p.value)

  write.table(cors, file=sprintf("autocor_%s.txt", transition.name), sep="\t", quote=F)
}



calc.autocors <- function(d, n.cores=1, B=100, lag=25)
{
  library(pbapply)
  if(n.cores==1)
  {
    rv = data.frame(do.call(rbind, pbapply(d, 1, autocor, B=B, lag=lag)))
    rv$ac = unlist(rv$ac)
    rv$p.value = unlist(rv$p.value)
  }else{
    #library(parallel)
    rv <- pblapply(1:nrow(d), function(i) autocor(as.numeric(d[i,]), B=B, lag=lag), cl = n.cores)
    rv = data.frame(do.call(rbind, rv))
    rv$ac = unlist(rv$ac)
    rv$p.value = unlist(rv$p.value)
    #registerDoSEQ()
  }
  rownames(rv) = rownames(d)
  return(rv)
}

#calculate auto-correlaton to find smoothly varying genes (probably need to generate a p-value)
autocor <- function(x, lag=1, B=0)
{
  ac = acf(x, plot = F, lag.max=lag+1)[[1]][lag] 
  if(B > 0)
  {
    check = 1:B
    for(i in 1:B)
    {
      check[i] = abs(acf(sample(x), plot = F, lag.max=lag+1)[[1]][lag]) > ac
      #cat(".")
    }
    p = sum(check)/B
  }else{
    p = NA
  }
  return(list("p.value"=p, "ac"=ac))
}



smooth_heatmap = function(obj, cells.use, pt, transition.name, genes)
{
  d = GetAssayData(obj)[genes, cells.use]
  du = smooth.data(du, loess.span=loess.span) 
  du = t(scale(t(du)))
  ## clip values 
  du[du > disp.range] = disp.range
  du[du< -disp.range] = -disp.range
  pseudo.cols = list("pseudo-time"=if(reverse) rev(pseudotime.cols) else pseudotime.cols)

  info(sprintf("Drawing transition [%s] heatmap for %s genes, %s cells", transition.name ,nrow(du), ncol(du)))
  write.table(rownames(du), file=sprintf("%s_genes_ordered.txt", transition.name), sep="\t", quote=F)
  aheatmap(du, filename=sprintf("%s_autocor.pdf", transition.name), main=sprintf("%s [%s cells]", transition.name, length(cells.use)),
    width=pdf.width, height=pdf.height, Colv = NA, Rowv = NA, annCol = ac, annColors = pseudo.cols, 
    color = colors, cexRow = cex.row, annRow = hb, cexCol = 0)
}






test_par = function(vec=1:100, ncores=4, progress=T)
{
    ### Start testing using future_lapply and progressor
    
    options(future.globals.maxSize = +Inf)
    plan(multisession, workers=ncores)
    #handlers(global = TRUE)
    #handlers("progress") #, "beepr")
    handlers(handler_progress(format="[:bar] :percent :eta :message"))
    
    if(!progress){rv = future.apply::future_lapply(vec, fit_lm)}else{
    start = Sys.time()
    rv = progressr::with_progress({
      p = progressr::progressor(along = vec)
      y = future.apply::future_lapply(vec, function(i, future.seed=123, ...) {
        p(sprintf("Processing gene=%g", i))
        fit_lm(x=NULL, y=NULL, nfit_pts=100)        
      })
    })
    }
    
    futile.logger::flog.info(sprintf("Elapsed: %s", Sys.time() - start))
    plan(sequential)
    rv
}





test_par_minimal = function(xs=1:1000, ncores=4)
{
    ### Start testing using future_lapply and progressor
    library(MASS)
    data(Boston)

    # function - calculate the mse from a model fit on bootstrapped samples from the Boston dataset
    model.mse <- function(x) {
        id <- sample(1:nrow(Boston), 200, replace = T)
        mod <- lm(medv ~ ., data = Boston[id,])
        mse <- mean((fitted.values(mod) - Boston$medv[id])^2)
        return(mse)
    }

    # initialising the list to pass to the apply functions
    x.list <- sapply(1:10000, list)

    # detect the number of cores
    ncores = detectCores(logical=F)

    print("Running sequential")
    p = system.time(a <- lapply(x.list, model.mse))
    print(p)

    clust <- makeSNOWCluster(ncores)
    print("Running parallel parLapply")
    p = system.time({
        clusterExport(clust, "Boston")
        a <- parLapply(clust, x.list, model.mse)})
    print(p)

    print("Running parallel future_lapply")
    p = system.time({
        plan(cluster, workers=clust)
        a <- future.apply::future_lapply(x.list, model.mse)})
    print(p)

    print("Running parallel pbl_apply")
    p = system.time({
        a <- pblapply(x.list, model.mse, cl=clust)})
    print(p)
}















