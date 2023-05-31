

library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(png)
library(RColorBrewer)
library(limma)
library(magrittr)
library(gridExtra)
library(knitr)

### Function to run pseudo bulk differential expression using EdgeR-LRT
### @param obj is a Seurat object
### @interp_max to reduce the sensitivity of this test to one or fewer super high expressing cells, we replace the top interp_max cells with the mean
### e.g. de = run_pseudobulk(hdm, compare = "control.hdm", compare_ref = "control.control")
run_pseudobulk = function(obj, compare, compare_ref, idents_celltype = "cell_labels", idents_compare = "condition", idents_sample = "mouse", mincells=10, mincounts=1, interp_max=1)
{
    md = obj@meta.data
    # Parameter checking
    if(!all(c(idents_compare, idents_celltype, idents_sample) %in% colnames(md))){stop(sprintf("'idents_compare, idents_sample, idents_celltype' args must all be columns in obj@meta.data: (%s)", paste(unique(colnames(md)), collapse = ", ")))}
    if(!compare %in% md[, idents_compare]){stop(sprintf("'compare' arg must be a level of 'idents_compare' arg (%s)", paste(unique(md[,idents_compare]), collapse = ", ")))}
    if(!compare_ref %in% md[, idents_compare]){stop(sprintf("'compare_ref' arg must be a level of 'idents_compare' arg (%s)", paste(unique(md[,idents_compare]), collapse = ", ")))}
    
    # Set up a singleCellExperiment object for EdgeR
    cat(sprintf("Running pseudobulk DE [%s vs %s]\n", compare, compare_ref))
    counts = GetAssayData(object = obj, slot = "counts", assay="RNA")
    md[, idents_compare] = make.names(md[, idents_compare]) #things with dashes like 'foxi1-ko' can mess up EdgeR
    compare = make.names(compare)
    compare_ref = make.names(compare_ref)
    sce_all = SingleCellExperiment(assays = list(counts = counts), colData = md) # Create single cell experiment object
    
    # Compute and add number of cells per sample
    sample_ids = as.factor(md[, idents_sample])
    n_cells = table(sample_ids) %>%  as.vector() ## Determine the 
    names(n_cells) = names(table(sample_ids))
    m = match(names(n_cells), sample_ids) ## Match the named vector with metadata to combine it
    cat(sprintf("Removing low-expressed genes, less than %s cells with at least %s counts.. \n", mincells, mincounts))
    
    # for each cell type, run edgeR w/ default parameters
    res = pbapply::pblapply(unique(md[, idents_celltype]), function(k) { 
        cat(sprintf("Processing %s.. ", k))
        sce = sce_all[, colData(sce_all)[, idents_celltype] == k]
        sce = sce[rowSums(counts(sce) > mincounts) >= mincells, ]
        cat(sprintf("%s genes tested. \n", nrow(sce)))
        ei = colData(sce) %>% data.frame() %>% dplyr::filter(!!as.symbol(idents_celltype) == k) # experiment info
        groups = colData(sce)[, c(idents_celltype, idents_sample)]  # Subset metadata to only include the cluster and sample IDs to aggregate across
        groups[, idents_sample] = factor(groups[, idents_sample])
        counts = counts(sce)
        if(interp_max > 0){
            counts = apply(counts, 1, function(x){x[order(x, decreasing = T)[1:interp_max]] = mean(x); x}) # Replace top interp_max values for each gene with the average. [FILTER STEP]'
        }else{counts = t(counts)}
        ac = aggregate.Matrix(counts, groupings = groups, fun = "sum")  #  Each row corresponds to aggregate counts for a cluster-sample combo
        if(all(make.names(c(compare, compare_ref)) %in% ei[, idents_compare])){
            m = match(names(n_cells), ei[, idents_sample]) 
            ei = data.frame(ei[m, ], n_cells, row.names = NULL) %>% dplyr::select(idents_sample, idents_compare, "n_cells")
            ei = na.omit(ei)
            design = model.matrix(~ 0 + ei[, idents_compare]) %>% 
              set_rownames(ei[, idents_sample]) %>% 
              set_colnames(levels(factor(ei[, idents_compare])))
            contrast = makeContrasts(contrasts = paste(compare, compare_ref, sep="-"), levels = design)
            celltype = str_split_fixed(rownames(ac), "_", 2)[,1]
            y = t(ac[celltype==k, ])
            y = DGEList(y, remove.zeros = TRUE)
            y = calcNormFactors(y)
            y = estimateDisp(y, design)
            fit = glmQLFit(y, design)
            fit = glmLRT(fit, contrast = contrast) # use the LRT test, it performed best in the 'confronting false discoveries' paper
            topTags(fit, n = Inf, sort.by = "none")$table %>% 
              dplyr::mutate(gene = rownames(.), celltype = k, contrast=paste(compare, compare_ref, sep="-")) %>% 
              dplyr::rename(p_val = PValue, p_adj = FDR)
        }else{
            cat(sprintf("WARN: Skipping %s, compare or reference groups are empty \n", k))
            NULL
        }
    })
    res = bind_rows(res)
    rownames(res) = NULL
    res %>% arrange(-abs(logFC))
}

