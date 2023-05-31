
#!/usr/bin/Rscript
###
# Single cell RNA-seq pipeline, May20
###
cat("Initialising\n")
library(scales)
library(methods)
library(getopt)
library(tools)
msg.trap <- capture.output( suppressMessages( source("util.r") ) ) # no warnings!

#require(DESeq2)


saved.1 = "pipeline.stage.1.Rdata"
saved.2 = "pipeline.stage.2.Rdata"
saved.3 = "pipeline.stage.3.Rdata"
saved.final = "pipeline.stage.final.Rdata"

show.settings <- function(opts)
{
    info("\n\n")
    info("----------- Single cell pipeline settings ----------")
    if(!is.null(opts$batch))
    {
        info(sprintf("   + Batch file                  %s",opts$batch))
    }else
    {
         info(sprintf("   + Counts table                 %s",opts$data))
         info(sprintf("   + QC metrics                  %s",opts$qc))
    }
    if(opt$ngenes > 0)
    {
        info(sprintf("   + Read only genes              %s",opts$ngenes))
    }
    if(opt$ncells > 0)
    {
        info(sprintf("   + Read only cells              %s",opts$ncells))
    }

    info(sprintf("   + Minimum genes per cell              %s", opts$min_genes))
    info(sprintf("   + Minimum cells per gene              %s", opts$min_cells))
    info(sprintf("   + Minimum transcriptome mapping      %s", opts$min_trans_mapped))
    if(opt$jackstraw)
    {
       info(sprintf("   + Run jackStraw test for sig. PCs  %s",opts$jackstraw))
    }
   
    info(sprintf("   + Find markers (SCDE)             %s",opts$markers))
    if(opt$markers)
    {
        info(sprintf("   + Compare to most sim cluster:   %s", opts$compare_sim_cluster))
        if(opt$use_pairwise_corr)
        {
            info(sprintf("   + Use pairwise corr to find most sim.: %s", opt$use_pairwise_corr))
        }
        if(opt$sim_cluster_min_size > 0)
        {
            info(sprintf("   + Most sim cluster has min size:        %s",opts$sim_cluster_min_size))
        }
    }
    info(sprintf("   + SCDE batch correct           %s", opts$scde_batch))
    info(sprintf("   + Many tSNE runs              %s", opts$many_tsne))
    if(opt$subsample > 0)
    {
        info(sprintf("   + Subsample                    %s",opts$subsample))
    }
    info(sprintf("   + Output dir                  %s",opts$output))
    if(!is.null(opt$drop_genes))
    {
        info(sprintf("   + Counts table                 %s",opts$drop_genes))
    }
    if(!is.null(opt$names_field_batch))
    {
        info(sprintf("   + Look for batch names in field %s",opts$names_field_batch))
    }
    if(!is.null(opt$tenx))
    {
        info(sprintf("   + 10X is %s ", opts$tenx))
    }
    if(!is.null(opt$tsne_precalc))
    {
        info(sprintf("   + Using tsne rotations         %s",opts$tsne_precalc))
    }else{
        info(sprintf("   + tSNE iterations              %s",opts$tsne_iters))
    }
    if(!is.null(opt$cluster_precalc))
    {
        info(sprintf("   + Using given clusters         %s",opts$cluster_precalc))
    }
    info(sprintf("   + P-value clustering            %s",opts$pvclust))
    info(sprintf("   + Analyse how many clusters     %s",opts$kclusters))
     if(!is.null(opt$cluster_interpret))
    {
        info(sprintf("   + Cluster interpretations        %s",opts$cluster_interpret))
    }
    info("----------------------------------------------------\n")
    info("\n\n")
}

#cat("Setting defaults\n")
#-------------- Default settings ----------------
do.find.bimodal.genes        = FALSE
use.hclust.assignments       = FALSE
rename_cells 				         = TRUE
run.seurat 					         = TRUE
run.seurat.cluster 	         = TRUE
load.complete.table          = TRUE # if this is true, extract the counts from a merged RSEM output table containing FPKM, TPM and expected_count, otherwise, 
                                    # load a counts matrix of and don't ask questions about the units.
pca.heatmaps                 = TRUE
n.pcs                        = 10
n.pc.extreme.genes           = 100
tsne.iterations              = 1000
min.genes                    = 1000
min.transcriptome.mapped     = 10
min.cells                    = 1
is.expr                      = 1 
do.plot                      = TRUE
dropped.samples              = "DROPPED" # samples in small clusters, that are dropped, are marked with this symbol
drop.small.clusters          = F
already.normalised           = FALSE
draw.variable.genes.heatmap  = TRUE
dseq                         = TRUE #if false, use the edgeR counts produced by the python pipeline.
use.tpm                      = T #if false use the expected_count
run.go.analysis              = FALSE
p.threshold.for.GO           = 0.0001
run.qc.analysis              = TRUE
run.jackstraw                = FALSE
dataset.name                 = NULL
visualise.coefficient.of.var = FALSE
annotate.marker.genes        = TRUE # add GO annotations to the markers for each cluster. is slow RN
trim.sample.names            = TRUE  # the qc data has no count string. so we trim it off, and tell seurat to look
                                    # for the group id in the first "_" delimited field, eg. Lgr5_A55_55 -> Lgr5

tsne.precalculated             = NULL  #can provide precalculated tSNE rotations for this dataset.
run.sparse.clustering          = FALSE
pipeline.many.tsne             = FALSE
pipeline.read.only.n.genes     = 0
pipeline.read.only.cells       = 0 
pipeline.scde.distance.metrics = TRUE
pipeline.remove.cell.cycle     = FALSE
pipeline.use.batch.file        = FALSE
pipeline.verbose               = FALSE
pipeline.build.scde.error.models = TRUE
pipeline.output.dir            = "SC_pipe"
find.markers.min.cluster.size  = 0 # if nonzero, small clusters (smaller than this param) are ignored when comparing to the most similar cluster.
find.markers.compare.sim.cluster   = FALSE # if this is false, markers are derived from only comparison to the background, if true, a comparison to the most similar cluster is also made
find.markers                   = FALSE # run SCDE to find genes enriched in clusters, takes ~ 2 hrs, should run on multi-core machine
find.markers.sample.n.cells    = 0 # if > 0, scde will subsample n cells from each group to compare. makes calculations much faster but doesn't consider all cells. debug only p. much

variable.genes.only            = FALSE # lost tlr4 once, probably not going to do this 
min.variance                   = 2  # keep only genes with variance > this val.
                                    # can also use the relative value of variance threshold to consider genes. i.e if its is 0.3 then the
                                    # bottom 30% of genes (ranked by variance across all samples) are dropped

# if(trim.sample.names)
# {
#     batch.names.field = 1
# }else
# {
#     batch.names.field = 3
# }

# cat("WARN: Setting names field to 3\n")
#------------------------------------------------


#cat("Setting command line args\n")

#----------- Command line args ------------------
#Column 1: the long flag name. A multi-character string.
# Column 2: short flag alias of Column 1. A single-character string.
# Column 3: Argument mask of the flag. An integer. Possible values: 0=no
# argument, 1=required argument, 2=optional argument.
# Column 4: Data type to which the flag’s argument shall be cast using storage.mode.
# A multi-character string. This only considered for same-row Column
# 3 values of 1,2. Possible values: logical, integer, double, complex, character. If
# numeric is encountered then it will be converted to double.
# Column 5 (optional): A brief description of the purpose of the option
args <- commandArgs(trailingOnly = TRUE)
n_args = length(args)
spec = matrix(c(
                'names_field_batch', 'bnf', 2, "integer", 
                'verbose', 'v', 2, "logical",
                'tenx', 'tx', 2, "logical",
                'help' , 'h', 0, "logical",
                'ngenes' , 'n', 2, "integer",    
                'ncells' , 'nc', 2, "integer",    
                'use_pairwise_corr', 'pwc', 2, 'logical',
                'compare_sim_cluster', 'csc', 2, "logical",
                'sim_cluster_min_size', 'ms', 2, "integer",
                'subsample', 'ss', 2, "integer",
                'batch', 'b', 2, "character",
                'cluster_interpret', 'ci', 2, "character",
                'data', 'd', 2, "character",
                'qc', 'q', 2, "character",
                'markers','m', 2, "logical",
                'jackstraw','js', 2, "logical",
                'run_id','o',2,"character",
                'remove_cell_cycle', 'rcc', 2, "logical",
                'scde_dist', 'dist', 2, "logical",
                'scde_batch', 'sb', 2, "logical",
                'drop_genes', 'dg', 2, "character",
                'tsne_iters', 'ti', 2, "integer",
                'min_cells', 'mc', 2, "integer", 
                'min_genes', 'mg', 2, "integer",
                'min_trans_mapped', 'mtm', 2, "integer",
                'pvclust', 'p', 2, "logical",
                'kclusters', 'k', 2, "logical",
                'many_tsne', 'mt', 2, "logical",
                'tsne_precalc', 'tp', 2, "character", 
                'cluster_precalc', 'cp', 2, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);
# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) | n_args==0) {
    cat(getopt(spec, usage=TRUE));
    q(status=1);
}
#set defaults 
if ( is.null(opt$names_field_batch ) ) { opt$names_field_batch = 1}
if ( is.null(opt$verbose ) ) { opt$verbose = FALSE }
if ( is.null(opt$ngenes ) ) { opt$ngenes = 0 }
if ( is.null(opt$ncells ) ) { opt$ncells = 0 }
if ( is.null(opt$use_pairwise_corr ) ) { opt$use_pairwise_corr = FALSE }
if ( is.null(opt$compare_sim_cluster ) ) { opt$compare_sim_cluster = FALSE }
if ( is.null(opt$sim_cluster_min_size ) ) { opt$sim_cluster_min_size = 0 }
if ( is.null(opt$subsample ) ) { opt$subsample = 0 }
if ( is.null(opt$markers ) ) { opt$markers = FALSE }
if ( is.null(opt$jackstraw ) ) { opt$jackstraw = FALSE }else{run.jackstraw=TRUE}
if ( is.null(opt$scde_dist ) ) { opt$scde_dist = FALSE }
if ( is.null(opt$remove_cell_cycle ) ) { opt$remove_cell_cycle = FALSE }
if ( is.null(opt$run_id)) {dataset.name=NULL} else{dataset.name=opt$run_id}
if ( is.null(opt$tsne_iters)) {opt$tsne_iters=tsne.iterations}
if ( is.null(opt$pvclust)) {opt$pvclust = FALSE}
if ( is.null(opt$kclusters)) {opt$kclusters = FALSE}
if ( is.null(opt$scde_batch)) {opt$scde_batch = FALSE}
if ( is.null(opt$tenx)) {opt$tenx = FALSE}
if ( is.null(opt$many_tsne)) {opt$many_tsne = FALSE}
if ( is.null(opt$min_genes)) {opt$min_genes=min.genes}
if ( is.null(opt$min_trans_mapped)) {opt$min_trans_mapped=min.transcriptome.mapped}else{
  min.transcriptome.mapped=opt$min_trans_mapped}
if ( is.null(opt$min_genes)) {opt$min_genes=min.genes}else{
  min.genes=opt$min_genes}
if ( is.null(opt$min_cells)) {opt$min_cells=min.cells}else{
  min.cells=opt$min_cells}

#cat("Processing args\n")
is.tenx                         = opt$tenx
batch.names.field               = opt$names_field_batch
cluster.interpret.file          = if(!is.null(opt$cluster_interpret)){file_path_as_absolute(opt$cluster_interpret)}else{NULL}
cluster.interpret               = NULL
pipeline.use.batch.file         = !is.null(opt$batch)
batch.file                      = if(!is.null(opt$batch)){file_path_as_absolute(opt$batch)}else{NULL}

counts.table                    = if(!is.null(opt$data)){file_path_as_absolute(opt$data)}else{NULL}
drop.genes.file                 = if(!is.null(opt$drop_genes)){file_path_as_absolute(opt$drop_genes)}else{NULL}
qc.file                         = if(!is.null(opt$qc)){file_path_as_absolute(opt$qc)}else{NULL}
pipeline.verbose                = opt$verbose
pipeline.remove.cell.cycle      = opt$remove_cell_cycle
pipeline.scde.distance.metrics  = opt$scde_dist
pipeline.output.dir             = sprintf("pipe_run_%s", dataset.name) 
pipeline.read.only.n.genes      = opt$ngenes
pipeline.read.only.n.cells      = opt$ncells
pipeline.pvclust                = opt$pvclust
pipeline.how.many.clusters      = opt$kclusters
pipeline.many.tsne              = opt$many_tsne
pipeline.find.markers           = opt$markers
find.markers.sample.n.cells     = opt$subsample
find.markers.min.cluster.size   = opt$sim_cluster_min_size
find.markers.compare.sim.cluster = opt$compare_sim_cluster
find.markers.use.pairwise.corr  = opt$use_pairwise_corr  # if this is false (default) the background markers are used to find the most similar cluster, if true pairwise corr is used
tsne.iterations                 = opt$tsne_iters
tsne.precalculated              = if(!is.null(opt$tsne_precalc)){file_path_as_absolute(opt$tsne_precalc)}else{NULL}
cluster.precalculated           = if(!is.null(opt$cluster_precalc)){file_path_as_absolute(opt$cluster_precalc)}else{NULL}
scde.batch.correct              = opt$scde_batch


if(is.null(dataset.name))
{
    error("No dataset name given!")
    stop("Provide a dataset name! --run_id")
}

if(!pipeline.use.batch.file)
{
    if(is.null(opt$data))
    {
        stop("Need either a counts table file (--data) or a batch file (--batch) arg!")
    }
    batch.groups = NULL
}
print(pipeline.output.dir)
#--------------------------------------------------------
info(sprintf("Creating output directory: %s", pipeline.output.dir))
dir.create(pipeline.output.dir, showWarnings = FALSE)
pipeline.output.dir = file_path_as_absolute(pipeline.output.dir)
setwd(pipeline.output.dir)
log.file = paste(getwd(), sprintf("pipeline_run_%s.log", dataset.name), sep="/")
log.to.console.and.file(log.file)
info(sprintf("Log file: %s", log.file))

# Other scripts:
info("Sourcing scripts")
scripts = c(
            "examine_subpopulation.r", 
      			"seurat.r",
			      "samples.r",
            "go_analysis.r",
      			"cluster.r",
            "pca.r",
            "qc.r",
            "tsne.r",
            "scde.r",
            "marker.r",
            "bimodal.r",
            "heatmap.r",
      			# "sclvm.r", # source this only if needed to avoid loading unnecessary libraries
            "graph.r"
            )
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
this.script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dir <- dirname(this.script.name)
initial_wd = getwd()


for (i in 1:length(scripts)) {
  	script.name = scripts[[i]]
  	other.name <- paste(sep="/", script.dir, script.name)
    print(paste("Sourcing", other.name, "from", this.script.name))
  	source(other.name)
}


info("Libraries loaded successfully.")
info("Starting pipeline..")
show.settings(opt)

info(sprintf("Dataset name: %s", dataset.name))

#Stage 1, load data 
#--------------------------------------------------------------------------
if(!file.exists(saved.1))
{
    if(pipeline.use.batch.file) # batch file specifies the locations and names of multiple datasets (and their QC files) to be merged
    {
          if(!is.null(cluster.interpret.file))
          {
              info(sprintf("Reading cluster interpretations from: %s", cluster.interpret.file))
              cluster.interpret=read.delim(cluster.interpret.file)
          }else
          {
              cluster.interpret = NULL
          }
          info("Batch file is ON..")
          if(is.tenx){use.tpm=T}
          batch.data = load.batch(batch.file, 
				                          use.tpm=use.tpm,
                                  parallel.read=T,
                                  min.genes=min.genes, 
                                  min.cells=min.cells, 
                                  min.trans.mapped=min.transcriptome.mapped,
                                  is.expr=is.expr,
                                  n.rows.only=pipeline.read.only.n.genes, 
                                  n.samples.only=pipeline.read.only.n.cells, plot.draw=T)   
          norm.counts = batch.data[["norm.counts"]]
          raw.counts = batch.data[["raw.counts"]]
          qc = batch.data[["qc"]]
          batch.groups = batch.data[["groups"]]
          complete_data = batch.data[["data"]]

          # eventually we need to do something clever about recognising that Lgr5Lo and Lgr5lo are the same
          # for now we convert sample names to upper case.
          info("Converting sample names to upper case ..")
          colnames(norm.counts) = toupper(colnames(norm.counts))
          colnames(raw.counts) = toupper(colnames(raw.counts))
          rownames(qc) = toupper(rownames(qc))
          qc$Sample = toupper(qc$Sample)
          
          info("Batch loading complete. Starting analysis..")

    }else{
          qc = NULL
    
          if(load.complete.table){


              if(use.tpm)
              {
                  info("Use TPM is ON")
                  count_string = "TPM"
              }else
              {
                  info("Use TPM is OFF (using normalised count)")
                  count_string = "expected_count" 
              }

              exclude = c("background","population", "nocell", "no.cell", "no_cell")
              
              # Load RSEM counts table
              expression_data = load_counts(counts.table, 
                                              count_string, 
                                              exclude, 
                                              use.row.names=T,
                                              n.rows.only=pipeline.read.only.n.genes,
                                              trim.sample.names=trim.sample.names)
              
              info("Counts loaded.")
              complete_data = expression_data[["data"]]
              

              # if(dseq){

              raw.counts = expression_data[["counts"]] # load counts data for the samples 
              info("Raw counts data size:")
        		  print(dim(raw.counts))
              if(is.tenx)
              {
                  warn("in 10X data mode")
                  
                  norm.counts = tpm(raw.counts)
              }else
              {
                  if(use.tpm)
                  {
                      info("Use TPM is ON ")
                      warn("DESeq normalisation is turned OFF. \'norm.counts\' is the same as raw.counts")
                      norm.counts = raw.counts
                  }else
                  {
                      info("Use TPM is OFF (Normalising)")
                      norm.counts = normalise.dseq(raw.counts)
                  }
              }
              
              extract.field=function(string,field=1,delim="_") return(strsplit(string,delim)[[1]][field])
              facs_type=factor(unlist(lapply(colnames(norm.counts), extract.field, batch.names.field,"_")))
              if(!is.null(qc.file)){
                  qc = load_qc(qc.file)
                  info(sprintf("Filtering. Keeping cells with more than %s genes detected above %s counts", min.genes, is.expr))
                  info(sprintf("    Keep genes expressed above %s counts in at least %s cells..", is.expr, min.cells))
                  norm.counts = filter.counts(norm.counts, 
                                              min.genes=min.genes, min.cells=min.cells, 
                                              is.expr=is.expr, draw.plot=T, qc=qc, colour.by=facs_type, 
                                              min.trans.mapped=min.transcriptome.mapped)
              }else{
                  qc = NULL
                  info(sprintf("Filtering. Keeping cells with more than %s genes detected above %s counts", min.genes, is.expr))
                  info(sprintf("    Keep genes expressed above %s counts in at least %s cells..", is.expr, min.cells))
                  print(corner(norm.counts))
                  if(min.cells > 0 | min.genes > 0)
                  {
                      norm.counts = filter.counts(norm.counts, min.genes=min.genes, 
                                              min.cells=min.cells, is.expr=is.expr, 
                                              colour.by=facs_type, draw.plot=T)
                  }else
                  {
                      info("min.cells and min.genes are 0. Not filtering..")
                  }
                  
              }
              raw.counts = raw.counts[rownames(norm.counts), colnames(norm.counts)]
                  
                  
              # }else{
              #     norm.counts = expression_data[["counts"]] # load counts data for the samples
                  
              #     if(!is.null(complete_data))
              #     {
              #         raw.counts = get_raw_counts(norm.counts, complete_data)

              #     }else{
              #         warn("Raw counts not provided!")
              #         raw.counts = NULL

              #         #print(head(norm.counts[, 1:20]))
              #     }
                  
              # }


          }else{
              
              info("Loading counts matrix ")
              gene.column = "GENE_SYMBOL"
              cat(sprintf("Looking for counts in %s.. \n", counts.table))
              cat(sprintf("Looking for gene names in column named: %s \n", gene.column))
              raw.counts = read.delim(counts.table)
              rownames(raw.counts) = raw.counts$GENE_SYMBOL
              raw.counts$GENE_SYMBOL <- NULL
              raw.counts = exp(raw.counts)

              if(already.normalised)
              {
                  info("Treating loaded counts as already normalised even though they are are called \'expected_count\'!")
                  norm.counts = raw.counts
              }
          }

          if(!is.null(drop.genes.file))
          {
              print("Dropping genes in file: ")
              print(drop.genes.file)
              print("Nrows before:")
              print(nrow(norm.counts))
              norm.counts = drop.genes.in.file(norm.counts, drop.genes.file)
              print("Nrows after:")
              print(nrow(norm.counts))
          }
    }

    if(variable.genes.only)
    {
        norm.counts= drop.non.varying.genes(norm.counts, min.var=min.variance)
        
    }
    info("Pipeline stage 1 complete. Saving..")
    save(list=c(
              "complete_data",
              "batch.groups",
              "norm.counts", 
              "raw.counts", "qc"), file=saved.1)

  }else{

    info("Found saved stage 1 object, loading instead of rerunning..")
    load(saved.1)
  }


  #---------- stage 2, initial analysis -----------------------
  if(!file.exists(saved.2)){

    #use SCDE to try to correct for the probability of dropouts
    if(pipeline.scde.distance.metrics)
    {
        # #debug:
        # keep.only.cells = 20
        # keep.only.genes = 5000
        # cat(sprintf("Only keeping %i cells ..\n",keep.only.cells))
        # raw.counts = raw.counts[, sample(ncol(raw.counts), keep.only.cells)]
        # if(keep.only.genes < nrow(raw.counts))
        # {
        #     raw.counts = raw.counts[sample(nrow(raw.counts), keep.only.genes), ]
        # }
        
        metrics = c("reciprocal", "mode.rel","direct")
        info("Running tSNE with SCDE-adjust metrics")
        print(metrics)
        initial_wd = getwd()
        dirname = "tSNE with SCDE-adjusted distance"
        dir.create(dirname, showWarnings = FALSE)
        setwd(dirname)
        scde.models = scde.compute.error.models(raw.counts)
        o.ifm = scde.models[["error.models"]]
        o.prior = scde.models[["expression.prior"]]

        plot.dropout.prob = T
        if(plot.dropout.prob)
        {
            scde.plot.dropout(scde.models[["error.models"]], scde.models[["expression.prior"]])
        }
        
        scde_wd = getwd()
        for(i in 1:length(metrics))
        {
            dist.type = metrics[i]
            dir.create(dist.type , showWarnings = FALSE)
            setwd(dist.type )
            
            dist = scde.dist(scde.models[["counts"]], 
                              scde.models[["error.models"]], 
                              scde.models[["expression.prior"]],
                              dist.type=dist.type)
            do.tsne(dist, n=20, parallel=T, max.iterations=tsne.iterations)
            setwd(scde_wd)
        }
        setwd(initial_wd)
    }else{
        o.ifm = NULL
        o.prior = NULL
    }

    


    if(pipeline.remove.cell.cycle)
    {
        source("sclvm.r")
        warn("Dropping genes that don't have an ensembl id..")
      	keep.genes = complete_data$ENSEMBL_ID != ""
      	warn(sprintf("%i genes left.. ", sum(keep.genes)))

      	# norm.counts = na.omit(norm.counts)
      	rownames(norm.counts) = complete_data$ENSEMBL_ID
      	norm.counts = norm.counts[keep.genes,]
      	raw.counts = raw.counts[keep.genes, ]
      	complete_data = complete_data[keep.genes, ]
      	
        info(sprintf("Ranges \n"))
        cat("Norm.Counts:\n")
        cat(range(norm.counts))
        cat("\n Raw.Counts:\n")
        cat(range(raw.counts))



        cat(sprintf("\n Correcting for cell-cycle using scLVM..\n"))
      	norm.counts = t(run_sclvm(norm.counts, 
                      							remove.cell.cycle=TRUE
                      							#draw.plots=FALSE
                                    ))
        cat("\n Corrected Norm.Counts: \n")
        cat(range(norm.counts))

      	# cat(sprintf("Size of returned (corrected) counts matrix: \n"))
      	# print(dim(norm.counts))
      	# print("Original size;")
      	# print(dim(raw.counts))
      	rownames(norm.counts) = rownames(raw.counts)
      	colnames(norm.counts) = colnames(raw.counts)
      	
      	# print(head(norm.counts))
        warn("Replacing negative values with 0's..")
        norm.counts[norm.counts < 0] <- 0

      	corr_counts_file = "corrected_counts.txt"
      	cat(sprintf("\n Writing corrected counts to %s ..\n", corr_counts_file))

      	write.table(norm.counts, file=corr_counts_file, sep="\t", quote=F)
      	#print correlation matrix after:
    }



    

    # load seurat:
    if(run.seurat){
        info("Size of counts dataframe:")
        info(paste(dim(norm.counts), collapse=", "))
        info("Starting seurat..")

        #------ Load seurat -----------------------------------------------#
        seurat.obj = seurat.setup(norm.counts, 
                                    pca.plots             = T,
                                    take.log              = T,
                                    run.jackstraw         = run.jackstraw, # examine significance of PCs, slow
                                    min.cells             = 1, #filtering already performed 
                                    min.genes             = 1,  #filtering arleady performed
                                    reduce.dims.for.tsne  = 6,
                                    use.ica.for.tsne      = FALSE,
                                    tsne.rotations        = tsne.precalculated, # this is NULL unless rotations are specified with --tsne_precalc
                                    tsne.iterations       = tsne.iterations,
                                    is.expr               = is.expr,
                                    seurat.names.field    = batch.names.field, 
                                    do.plot               = do.plot)
        #------------------------------------------------------------------#

        info("Seurat cell names:")
        info(head(seurat.obj@cell.names))

        info("Adding raw counts to seurat.obj:")
        seurat.obj@raw.data = data.frame(raw.counts)

        info("Building cell ID matrix:")
        cell.types = as.character(unlist(fetch.data(seurat.obj, c("orig.ident"))))

        if(!is.null(batch.groups))
        {
            ident = cbind(cell.types, as.character(batch.groups))
            colnames(ident) = c("Cell Type", "Batch")
            rownames(ident) = seurat.obj@cell.names
        }else
        {
            info("Building ID matrix..")
            ident = matrix(cell.types)
            rownames(ident) = seurat.obj@cell.names
            colnames(ident) = "Cell Type"
        }
        info("ID types:")
        info(paste("       ", colnames(ident)))
        seurat.obj = seurat.cluster(seurat.obj, find.markers=T, dbclust.eps = 30)
        
        info("Seurat clustering complete.") 
        sigs = paste(path.expand("~"),"/ref/gene_lists/gut_circuits/signatures/Gene.Lists.all.updated.June2016.txt", sep="")
        info("Scoring cells based on known gene signatures..")
        info(sprintf("Signature file: %s", sigs))
        gene.lists = read.delim(sigs)
        scores = score.cells(seurat.obj@data, gene.lists)
        seurat.obj = set.attribute(seurat.obj, scores)
	      seurat.obj@data.info[,"Batch"] = batch.groups
        
        # if(draw.variable.genes.heatmap){
            
        #     # if qc is provided, this should call Sam R's heatmap code, and 
        #     # show the expression of the most variable genes and also the
        #     # most relevant QC metric (number of genes)
        #     info("Drawing heatmap of expression of variable genes detected by Seurat..")
        #     setwd(pipeline.output.dir)
           
        #     d = heatmap.signature(seurat.obj, var.only=T, pdf.output="variable_genes_heatmap.pdf")
           
        #     pdf("variable_genes_heatmap_unscaled.pdf", width=10, height=12)
        #     d = heatmap.signature(seurat.obj, var.only=T, scale=F)
            
        #     counts.var = d$data 
        #     dev.off()
            
        #     exp.fit.variable.genes = 100
        #     pdf("variable_genes_expfit.pdf")
        #     var.genes.exp.fit = get.variable.genes(norm.counts, show.n = exp.fit.variable.genes, winsorize=T, fit.all=T)
        #     dev.off()
        #     pdf("var.genes.expfit.heatmap.pdf", width=12, height=18);
        #     aheatmap(norm.counts[head(var.genes.exp.fit$Gene, n=exp.fit.variable.genes),]) 
        #     dev.off()
        #     info("Variable genes heatmap complete.")
        # }

       

        #print(seurat.obj)
        # only keep the cells and genes that pass Seurats filter, this
        # makes sure that all the matrices are the right size etc.

        # this means norm.counts are scaled and normalised and log transformed
        # from here on in.
        norm.counts = seurat.obj@data

         # run sparse clustering:
        if(run.sparse.clustering)
        {
            m = seurat.obj@scale.data[seurat.obj@var.genes, ]
            shc = sparse.hclust(t(m))
            sparse.genes = rownames(m)[which(shc$ws > 0)]
            info(sprintf("Hclust finished in sparse basis of %s genes. ", length(sparse.genes)))
            aheatmap(norm.counts[sparse.genes,], filename="sparse.genes.pdf")
        }
        

    }else
    {
        seurat.obj = NULL
        norm.counts = log2(norm.counts + 1)
    }

    if(pipeline.many.tsne)
    {
        if(tsne.iterations < 2000)
        {
            warn("Increasing tSNE iterations to 2000")
            tsne.iterations = 2000
        }
        many.tsne(seurat.obj, 
              			iters=tsne.iterations, 
              			n.var.genes=2000, 
              			perp.values=c(8, 10, 12, 15, 20, 25), 
              			pc.values=c(5:15), 
              			n.runs=12, 
              			barnes.hut=T, 
              			n.cores=4, 
              			snow=F)
    }


    if(is.null(raw.counts) & !is.null(complete_data))
    {
        raw.counts = get_raw_counts(norm.counts, complete_data)
    }

    setwd(pipeline.output.dir)
    info("Pipeline stage 2 complete. Saving..")
    save(list=c(
              "o.ifm",
              "o.prior",
              "complete_data",
              "seurat.obj", 
              "ident", 
              "batch.groups",
              "norm.counts", 
              "raw.counts", "qc"), file=saved.2)
}else
{
    info("Found saved stage 2 object, loading instead of rerunning..")
    load(saved.2)
}

# Stage 2 complete. Save Seurat object, and if rerun, reload from here and 
# begin again.


scaled.log.counts= t(scale(t(log2(norm.counts+1)),center=TRUE,scale=TRUE))
if(pipeline.pvclust)
{
    info("Running P-value clustering..")
    pv = cluster.pvclust(scaled.log.counts, min.p = 0.05, nboot=1000, parallel = T, pdf.output = T)
    info("P-value clustering done.")
}

if(pipeline.how.many.clusters)
{
    info("Estimating the number of clusters in the data..")
    cluster.how.many(scaled.log.counts)
    info("Cluster size estimation done.")
}

if(!file.exists(saved.3))
{
    if(!is.null(cluster.interpret.file) & is.null(cluster.interpret))
    {
        info(sprintf("Reading cluster interpretations from: %s", cluster.interpret.file))
        cluster.interpret=read.delim(cluster.interpret.file)
    }else
    {
        if(is.null(cluster.interpret.file))
        {
            cluster.interpret= NULL
        }
    }
    info("Running hierarchical clustering (hclust)..")
    hclust.assignments = cluster.hierarchical(norm.counts, 
            take.log           = FALSE, # careful not to take log twice (seurat already did!)
            # use.manual.labels  = c("Lgr5+ [G1&S]", "Lgr5+ [G2&M]", 
            #                         "Paneth", "Enterocyte", "Enteroendocrine","Goblet","Stem", "Paneth-Goblet"),
            sample.ident       = ident,
            show.gene.type     = F,
            show.sample.type   = T,
            show.cluster.assignment = T,
            row.dendro         = F, 
            row.labels         = F, 
            min.cluster.size   = 15,
            heatmap.draw       = F,
            # heatmap.name       = "all_genes_heatmap.pdf",
            # heatmap.title      = paste("All genes \n", dataset.name),
            pdf.output         = T,
            verbose            = pipeline.verbose,
            raw.counts         = raw.counts)$hc.col


    
    #-------------- Choose the clusters to use ---------------
    if(!is.null(cluster.precalculated))
    {
        info(sprintf("Using clusters from %s ", cluster.precalculated))
        read.in = read.delim(cluster.precalculated)
        clusters = unlist(read.in)
        names(clusters) = rownames(read.in)
    }else
    {
        if(use.hclust.assignments)
        {
            info("Using hclust assignments!")
            clusters = hclust.assignments
        }else{
            info("Using tSNE/DBSCAN assignments!")
            clusters = seurat.obj@ident
        }
    }

    info("Running pamk..")
    pamk.best <- pamk(t(seurat.obj@data))
    cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n")

    info("Building PC distance matrix..")
    pcd = pc.dist(t(seurat.obj@data), pcs.use=1:10)

    info("Running graph clustering..")
    g = graph.cluster(t(seurat.obj@data), k=20, dm=pcd)

    seurat.obj@data.info$hclust.ID = factor(hclust.assignments)
    seurat.obj@data.info$tSNE_dbscan.ID = seurat.obj@ident
    seurat.obj@data.info$PAM.ID = factor(pamk.best$pamobject$clustering)
    seurat.obj@data.info$graph.cluster.ID = factor(g$partition)

    # print(clusters)
    # print(is.na(as.numeric(clusters)))
    if(!is.null(cluster.interpret))
    {
        info("Current clusters: ")
        info(table(clusters))
        info("Renaming clusters using interpretation file..")
        info(cluster.interpret)
        clusters = cluster.interpret$Interpretation[match(clusters, cluster.interpret$Cluster.ID)]

    }else
    {
        warn("No given cluster interpretation!")
    }
    clusters.names = unique(clusters)
    names(clusters) = seurat.obj@cell.names
    seurat.obj@data.info$DBclust.ident = clusters # so feature.plot.scale looks good. 
    
    info("Drawing correlation matrix..")
    cor.matrix(seurat.obj, title=paste(dataset.name, "correlation matrix"), clusters=clusters)
    #---------------------------------------------------------


    if(drop.small.clusters)
    {
        clusters = drop.small(norm.counts, 
                              clusters, 
                              dropped.samples=dropped.samples)
        clusters.names = unique(clusters)
        info("Size of counts dataframe:")
        info(paste(dim(norm.counts), collapse=", "))
    }else
    {
        info("Drop small clusters is OFF..")
    }

    if(visualise.coefficient.of.var)
    {
        vis.coeff.hists(log2(norm.counts+1), clusters)
    }

    if(do.find.bimodal.genes)
    {
        info("Finding bimodal genes in each cluster..")
        find.bimodal.genes(norm.counts, clusters, write.tables=T)
    }


    # info("Generating PCA bi-plots..")
    # pca = prcomp(t(norm.counts), scale=T, center=T)

    if(is.null(seurat.obj@pca)){
        info("Visualising PCA using groups:")
        groups = unlist(fetch.data(seurat.obj, c("orig.ident")))
        names(groups) = seurat.obj@cell.names
        info(table(groups))
        visualise.pca(seurat.obj@pca.obj, groups, n_components=6)
    }
    

    #Run QC analysis:
    if(run.qc.analysis)
    {
        if(is.null(qc))
        {
            if(!is.null(qc.file)){
                qc = load_qc(qc.file)
            }else{
                qc = NULL
            }
        }
        
        if(!is.null(qc))
        {
            setwd(pipeline.output.dir)
            info("QC dataframe size:")
            info(paste(dim(qc), collapse=", "))
            info("Keeping QC only for the samples we have counts for..")
            
            use.row.name = FALSE
            if(use.row.name)
            {
                info("Using QC rownames:")
                info(head(rownames(qc)))
                qc.cells = rownames(qc) %in% seurat.obj@cell.names
                qc = qc[qc.cells, ]
                qc$Sample = rownames(qc)
                info("QC dataframe size:")
                info(paste(dim(qc), collapse=", "))
                counts.not.qc = seurat.obj@cell.names[!(seurat.obj@cell.names %in% rownames(qc))]
            }else
            {
                info("Using QC sample names:")
                info(head(qc$Sample))
                qc.cells = qc$Sample %in% seurat.obj@cell.names
                qc = qc[qc.cells, ]
                info("QC dataframe size:")
                info(paste(dim(qc), collapse=", "))
                counts.not.qc = seurat.obj@cell.names[!(seurat.obj@cell.names %in% qc$Sample)]
            }
            
            info("These cells are in seurat.obj@cell.names but not QC:")
            info(counts.not.qc)

            
            check.pcs(qc, pca, groups)
            if(length(unique(clusters)) > 2)
            {
                check.qc.in.clusters(qc, clusters)
            }else
            {
                warn("No clusters detected. Cannot check QC for them")
            }
            
            info("Adding QC to seurat object..")
            seurat.obj = set.attribute(seurat.obj, qc)
            info("QC analysis complete.")
        }else{
            error("Failed to load QC!")
            check.pcs(NULL, pca, groups, counts=seurat.obj@data)
        }
    }


    if(!is.null(qc))
    {
        info("Embedding QC into seurat object..")
        seurat.obj = set.attribute(seurat.obj, qc)
        # note, seurat's 'ngene', in seurat.obj@data.info is the number of rows > is.expr 
        cluster.complexity(seurat.obj, clusters=clusters, pdf.output=T)
        reads.mapping.bar(seurat.obj)
    }
    

    if(pca.heatmaps)
    {
        info("Generating PC-Extreme Gene heatmaps..")
        info(sprintf("PCS:   %i", n.pcs))
        info(sprintf("Genes: %i", n.pc.extreme.genes))

        pc.heatmaps.dir = paste(pipeline.output.dir, "pc_extreme_heatmaps",sep="/")
        dir.create(pc.heatmaps.dir, showWarnings = FALSE)
        setwd(pc.heatmaps.dir)
      
        # info("Running PCA..")
        #pca = prcomp(t(norm.counts), scale=T, center=T)
        
        # rahul's got his labels backward
        loadings <- seurat.obj@pca.x #pca$rotation


        for(pc in 1: n.pcs){
            info(sprintf("Writing top loaded genes for PC-%i ..", pc))
            l = abs(loadings[, pc])
            top.genes = rownames(data.frame(head(sort(l, decreasing=TRUE), n=n.pc.extreme.genes)))
            #colnames(top_loadings) = c(paste("PC-",pc, sep=""))
            info(top.genes)
            info(sprintf("Drawing PC-%i heatmap ..",pc))
            aheatmap(seurat.obj@data[top.genes, ], scale="row", filename=paste("PC",pc,"hmap.pdf",sep="_"), color=c("Purple","Black", "yellow"))
            
            pdf(paste("PC",pc,"hmap.seurat.pdf",sep="_"))
            pcHeatmap(seurat.obj,pc.use = pc,do.balanced = FALSE)
            dev.off()
            # variable.clusters = cluster.hierarchical(seurat.obj@data, 
            #         genes          = top.genes, 
            #         heatmap.draw   = T,
            #         heatmap.name   = paste("PC",pc,"hmap.pdf",sep="_"),
            #         heatmap.title  = paste(sprintf("PC-%i genes", pc), dataset.name),
            #         color.palette.heatmap = c("Purple","Black", "yellow"), #RdBu", #"SamR",
            #         show.gene.type = F, 
            #         sample.ident   = ident,
            #         row.label.size = 0.4, 
            #         row.dendro     = F, 
            #         verbose        = pipeline.verbose,
            #         pdf.output     = T)         
        }
        setwd(pipeline.output.dir)

    }

    # set the cluster assignments (possibly generated by hclust) on the seurat object, so that
    # figures (feature and tsne plots) generated by seurat refer to the right 
    # clusters.
    if(!is.null(seurat.obj)){
        info("Labeling seurat objects with hclust labels..")
        seurat.obj = set.cluster.labels(seurat.obj, hclust.assignments)
	if(drop.small.clusters){
            info("Removing any cells in dropped clusters..")
            not.dropped = which(clusters != dropped.samples)
            clusters = clusters[not.dropped]
            if(!is.null(qc))
            {
                qc = qc[not.dropped, ]
            }
            norm.counts = norm.counts[ , not.dropped]
            raw.counts = raw.counts[ , not.dropped]
            seurat.obj = subsetData(seurat.obj, cells.use=colnames(norm.counts))
            # print("New seurat object (subsetted) has ident:")
            # print(length(seurat.obj@ident))
            # print(seurat.obj@ident)
        }

    }

    if(!is.null(seurat.obj)){
        if(do.plot)
        {
            pdf(paste("tsne_labeled_with_hclust.pdf"),width=11, height=8.5)
        }
        tsne.plot(seurat.obj, do.label = TRUE, label.pt.size = 1, pt.size=1)
        if(do.plot)
        {
            dev.off()
        }
  

        if(!use.hclust.assignments)
        {
            seurat.obj = set.cluster.labels(seurat.obj, clusters)
        }
        if(do.plot)
        {
            pdf(paste("tsne_labeled_with_dbscan.pdf"),width=11, height=8.5)
        }
        tsne.plot(seurat.obj, do.label = TRUE, label.pt.size = 1, pt.size=1)
        if(do.plot)
        {
            dev.off()
        }
    }


    info("Downsampling genes using hclust..")
    downsampled.counts = downsample.genes(norm.counts, n_genes = 500)
    pdf("downsampled.500.pdf", width = 16, height = 20);aheatmap(downsampled.counts, info = T, cexRow = 12, scale = "row"); dev.off()


    info("Pipeline stage 3 complete. Saving..")
    save(list=c(
              "o.ifm",
              "o.prior",
              "downsampled.counts",
              "complete_data",
              "seurat.obj", 
              "ident", 
              "clusters",
              "batch.groups",
              "norm.counts", 
              "raw.counts", "qc"), file=saved.3)
}else{
    info("Found saved stage 3 object, loading instead of rerunning..")
    load(saved.3)
}



# ------------------- Search for markers ----------------------
if(pipeline.find.markers){ 
    info("Clusters:")
    info(table(clusters))
    info("Searching for marker genes in each cluster:")
    if(find.markers.sample.n.cells > 0)
    {
        info("Subsampling is on.")
        info(sprintf("Will sample %i cells from each cluster..", find.markers.sample.n.cells))
    }

    
    clusters = check.cluster.names(clusters)  
    checkfile = "markers.ok"
    markers.dir = paste(pipeline.output.dir, "markers",sep="/")
    dir.create(markers.dir, showWarnings = FALSE)
    setwd(markers.dir)
    for(i in 1: length(unique(clusters)))
    {
        cluster.id = as.character(unique(clusters)[i])
        if(!identical(cluster.id, dropped.samples))
        {
            cat("\n \n")
            info(sprintf("Examining \'%s\', cluster %i of %i..", cluster.id, i, length(unique(clusters))))
            info(sprintf("Creating \'%s\' directory.. ", cluster.id))
            dir.create(cluster.id, showWarnings = FALSE)
            setwd(cluster.id)
            info(sprintf("Moving to %s", getwd()))
            if(!file.exists(checkfile))
            {
                if(scde.batch.correct)
                {
                    use.batch = batch.groups
                }else
                {
                    use.batch= NULL
                }
                # call SCDE to run differential expression tests on the cells in this cluster.
                cluster_info = examine_cluster(norm.counts, 
                                 #complete_data, 
                                 clusters, 
                                 use.pairwise.corr      = find.markers.use.pairwise.corr,
                                 most.similar.min.size  = find.markers.min.cluster.size,
                                 compare.most.similar   = find.markers.compare.sim.cluster,
                                 cluster.id             = cluster.id,
                                 annotate.marker.genes  = annotate.marker.genes,
                                 raw.counts             = raw.counts,
                                 sample.n.cells         = find.markers.sample.n.cells,
                                 batch                  = use.batch, 
                                 verbose                = pipeline.verbose,
                                 single_cell_data       = TRUE,
                                 seurat.obj             = seurat.obj)

                # draw an average expression heatmap of top n markers
                n.markers = 20
                if(!is.null(cluster_info))
                {
                    if(nrow(cluster_info) > 1)
                    {
                        plot.markers(seurat.obj=seurat.obj, de.sig=cluster_info, clusters=clusters)
                        if(nrow(cluster_info) < n.markers)
                        {
                            n.markers = nrow(cluster_info)
                        }
                        av = average.heatmap(norm.counts, 
                                      clusters, genes=unique(head(cluster_info$GENE_SYMBOL, n=n.markers)), scale=T,  
                                      pdf.output=TRUE, pdf.name=paste(cluster.id,".Top", n.markers,"Markers.Z.Score.pdf", sep=""))
                        av = average.heatmap(norm.counts, 
                                      clusters, genes=unique(head(cluster_info$GENE_SYMBOL, n=n.markers)), scale=F,  
                                      pdf.output=TRUE, pdf.name=paste(cluster.id,".Top", n.markers,"Markers.pdf", sep=""))
                    }else
                    {
                      warn(sprintf("No markers detected for %s", cluster.id))
                    }
                  }else
                  {
                      warn(sprintf("No markers detected for %s", cluster.id))
                  }
                  info(sprintf("Find markers complete for %s. ", cluster.id))
                  fileConn<-file(checkfile)
                  writeLines(c(cluster.id,"Find markers completed."), fileConn)
                  close(fileConn)
            }else
            {
                info("Markers already found for this cluster. Moving on..")
            }
            setwd(markers.dir)
        }else
        {
            cat(sprintf("Skipping %s !!\n", dirname))
        }
    }

    # add the average expression of each marker gene (in each cluster)
    add.average.expression(norm.counts, clusters)
    
    # #(re) draw feature plots of top_n markers
    # draw.feature.plots(norm.counts, clusters, top_n)

    #do markers heatmap:
    info("Drawing markers heatmaps..")
    draw.markers.heatmap(norm.counts,
                          cluster.names = unique(clusters), 
                          ident, 
			  cluster.assignments = clusters,
			  output.pdf = T)

    info(sprintf("Running GO analysis of clusters [p=%s]", p.threshold.for.GO))
    run.go.analysis(norm.counts, clusters, p.threshold=p.threshold.for.GO)
    

}else
{
    info("Examining clusters is OFF..")
}

#useful plot of the cluster labels of bulk (unsorted cells)
bulk.cell.ids(seurat.obj, clusters)

if(pipeline.build.scde.error.models & is.null(o.ifm) | is.null(o.prior))
{
    info(sprintf("Will fit SCDE models for %s cells", ncol(raw.counts)))
    #info("Not using groups to compute error models as per the bug Carl found")
    scde.models = scde.compute.error.models(raw.counts[, use.cells], verbose=verbose, n.cores=n.threads) #, groups=groups)
    o.ifm = scde.models[["error.models"]]
    counts_data = scde.models[["counts"]]
    o.prior = scde.models[["expression.prior"]]
    info(sprintf("Built error models for %s cells", ncol(counts_data)))
}

info("Pipeline complete. Saving..")
    save(list=c(
              "o.ifm",
              "o.prior",
              "downsampled.counts",
              "complete_data",
              "seurat.obj", 
              "ident", 
              "clusters",
              "norm.counts", 
              "raw.counts", "qc"), file=saved.final)
info("Exiting..")
# --------------------------------------------------------------





