
### Utilities functions to analyze bulk RNA-seq data with DESeq2
library(DESeq2)
library(EDASeq)
library(pheatmap)

setupSampleAnnotations <- function(cts, groups)
{
    info("Making colData")
    colData = data.frame(colnames(cts), check.names=F, stringsAsFactors=F)
    colnames(colData) = "Sample"
    colData$Condition = factor(groups)
    rownames(colData) = colData$Sample
    colData$Sample = NULL
    colData
}

setupDEseq2 <- function(counts, colData, min.reads=1, quantile.norm=T, design.formula=~Condition)
{

  info("Setting up DESeq2 object")
  
  if(nrow(colData) != ncol(counts)){stop("colData and counts don't have the right dimensions!")}
  info(sprintf("Data Dimensions: Counts -> %s, colData -> %s", 
    paste(dim(counts), collapse="x"), paste(dim(colData), collapse="x")))
  rounded = round(counts)
  rounded = apply(counts, 2, function(x) {storage.mode(x) <- 'integer'; x})
  
  print(colnames(rounded))

  # info("ColData:")
  # print(colData)
  # n.conds = length(levels(colData$Condition))
  # conds = paste(levels(colData$Condition), collapse=",")
  # info(sprintf("Setting up DESeq2 object with %s conditions: %s", n.conds, conds))
  # info(sprintf("Quantile normalisation is %s", if(quantile.norm){"ON"}else{"OFF"}))

  rownames(colData) = colnames(rounded)

  dds = DESeqDataSetFromMatrix(countData = rounded, colData = colData, design = design.formula) # ~Batch + Condition)
  dds <- dds[ rowSums(counts(dds)) > min.reads, ]

  info(sprintf("Sample names:"))
  print(colnames(counts(dds)))
  if(quantile.norm)
  {
    info("Using EDA quantile normalisation")
    eda.offset <- betweenLaneNormalization(as.matrix(round(counts(dds))),which="full",offset=TRUE)
    nf <- exp(-1*eda.offset)
    # print(dim(nf))
    # print(corner(nf))
    #nf <- nf[,c("ex1T6","EX2T7","ex3T7","ex1T2","EX2.T4T5","EX3.T4T5")];
    nf <- nf/mean(nf)
    #info("EDA norm factors:")
    # print(head(nf))
    # print(dim(nf))
    normalizationFactors(dds) <- nf;
    
  }else{
    info("Using DESeq2 normalisation")  
    #dds <- estimateSizeFactors(dds) # redundant, i think?
  }
  info("Running DEseq tests")
  #dds <- estimateDispersions(dds) # redundant, i think?
  #dds <- nbinomLRT(dds, reduced = ~ Batch)
  dds <- DESeq(dds, quiet=F) #, reduced = ~ exp)  
  return (dds)
}

ma.plot <- function(x, plot.name, p.thresh=0.1, label.size=5, y.axis.label=NULL, x.axis.label=NULL, cex.all=20, max.x=18)
{
  library(ggrepel)

  if(is.null(y.axis.label)){y.axis.label = "Log2 fold change (DESeq2)"}
  if(is.null(x.axis.label)){x.axis.label = "Mean expression\n(Log2 normalized transcript count, DESeq2)"}

  ggsave(ggplot(x, aes(x=log2(baseMean+1), y=log2FoldChange, label=gene)) + geom_point(alpha=0.2, size=0.2)  + theme_bw() + xlim(c(0, max.x)) + #ylim(c(-12,12)) + 
      geom_point(data=subset(x, padj<p.thresh), aes(colour=cell.type), alpha=1, colour="red") +
      #geom_text_repel(data=subset(x, padj<p.thresh), size = label.size, box.padding = unit(0.35, "lines"),
          #point.padding = unit(0.3, "lines"), segment.size = 0.1, max.iter=10e3, fontface = 'italic') + 
        ylab(y.axis.label) + xlab(x.axis.label) + geom_text(data=subset(x, padj<p.thresh), size = label.size, fontface = 'italic', vjust=2) +
      theme(text = element_text(size=cex.all, colour="gray22")), filename=sprintf("%s_MA.pdf", plot.name), width=10, height=9) 

  ggsave(ggplot(x, aes(x=log2(baseMean+1), y=log2FoldChange, label=gene)) + geom_point(alpha=0.2, size=0.2)  + theme_bw() + xlim(c(0, max.x)) + #ylim(c(-12,12)) + 
      geom_point(data=subset(x, padj<p.thresh), aes(size=-log10(padj)), alpha=0.4, colour="red") +
      #geom_text_repel(data=subset(x, padj<p.thresh), size = label.size, box.padding = unit(0.35, "lines"),
          #point.padding = unit(0.3, "lines"), segment.size = 0.1, max.iter=10e3, fontface = 'italic')+ 
          ylab(y.axis.label) + xlab(x.axis.label) + geom_text(data=subset(x, padj<p.thresh), size = label.size, fontface = 'italic', vjust=2) +
      theme(text = element_text(size=cex.all, colour="gray22")), filename=sprintf("%s_MA_pval_sized.pdf", plot.name), width=10, height=8) 
}


cdf.plot <- function(counts, cell.type, suffix="")
{
  # draw cumulative distribution factor to visualise sample-sample differences
  cdf = sprintf("%s_CDF_%s.pdf", cell.type, suffix)
  info(sprintf("Drawing CDF to %s", cdf))
  ctm = melt(counts)
  print(colnames(ctm))
  n.col = length(unique(ctm$variable))
  ggsave(ggplot(ctm,aes(log2(value+1),colour=variable)) + stat_ecdf(size=1) +  
    scale_colour_manual("Sample", values=default.cols(n.col)) + 
    ggtitle(sprintf("Cumulative distribution function [%s]", cell.type)) + theme_bw() +
    xlab("Log2 RSEM expected count") + ylab("") + 
    theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) + 
    theme(plot.title=element_text(vjust=-3)) + 
    theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        text = element_text(size=20, colour="gray22")) , filename=cdf, width=9, height=8)
}


### a wrapper around calling 'results' on a DESeq2 object that also draws some plots.
####
process.results <- function(dds, 
  p.thresh=0.05, 
  gene.enrich=F, 
  clip=2, 
  p.thresh.label=0.01, 
  label.size=5, 
  flip.direction=F, 
  look_for_extras=F,
  contrast=NULL)
{

  if(is.null(contrast))
  {
    warn("No contrast given. DESeq2 Manual: With no arguments
to results, the results will be for the last variable in the design formula, and if
this is a factor, the comparison will be the last level of this variable over the first
level.")
    variable = colnames(colData)[1]
    compare = rev(levels(colData)[,1])
  }else{
    if(length(contrast)!= 3){
      stop("DESeq2 Manual: constrast argument must be a character vector of length three: the name of the variable, the name of the factor
level for the numerator of the log2 ratio, and the name of the factor level for the denominator.")
    }
    variable = contrast[1]
    compare = contrast[2:3]
  }

  plot.name = paste(variable, compare[1], "vs", compare[2], sep="_")
  info(sprintf("Variable: %s [%s vs %s]", variable, compare[1], compare[2]))
  info("Extracting DE genes")
  res <- results(dds, contrast=c(variable, compare))
  
  x = data.frame(res)
  x$gene = rownames(x)
  
  
  sig = which(x$padj < p.thresh)
  sig.genes = rownames(x)[sig]
  n.sig = sum(x$padj < p.thresh, na.rm=TRUE)
  info(sprintf("%s genes significantly DE (adjusted p-val < %s)", n.sig, p.thresh))
  x <- x[order(x$padj),]

  info("Drawing MA plot")
  ma.plot(x, plot.name, p.thresh=p.thresh.label, label.size=label.size, y.axis.label=sprintf("Log2 fold change (%s DESeq2)", plot.name))
  print(summary(x))
  info("Top sig. genes:")
  
  print(head(x[, c("padj", "pvalue", "log2FoldChange")], n=100))


  if(flip.direction){x$log2FoldChange = -x$log2FoldChange}

  y = x[sig.genes, c("padj", "pvalue", "log2FoldChange")]
  y$Gene = rownames(y)
  y = y[order(abs(y$log2FoldChange), decreasing=T),]

  up = y[y$log2FoldChange>0,]
  down = y[y$log2FoldChange<0,]
  write.table(up, file=sprintf("%s_UpReg_FDR_%s.txt", plot.name, p.thresh), quote=F, sep="\t", row.names=F)
  write.table(down, file=sprintf("%s_DownReg_FDR_%s.txt", plot.name, p.thresh), quote=F, sep="\t", row.names=F)
  write.table(x[order(abs(x$padj)),], file=sprintf("%s_Complete.txt", plot.name), quote=F, sep="\t", row.names=F)

  # if(gene.enrich)
  # {
  #   background = as.character(rownames(norm.counts))
  #   go.analysis(as.character(up$Gene), background, sprintf("%s_Upreg_Sig_p_%s", plot.name, p.thresh), annotation="mm10")
  #   go.analysis(as.character(down$Gene), background, sprintf("%s_Downreg_Sig_p_%s", plot.name, p.thresh), annotation="mm10")
  # }

  norm.counts = counts(dds, normalized=T)
  data = log2(norm.counts[sig.genes,]+1)
  me = data.frame(Mean_Log2_TPM=rowMeans(data))
  show.data = t(scale(t(data)))
  show.data[show.data > clip] <- clip
  show.data[show.data < -clip] <- -clip


  cond = dds@colData[, variable]
  cond = data.frame(cond)
  colnames(cond) = "cond"
  pheatmap(show.data, color=colorRampPalette(rev(brewer.pal(11, "RdBu")))(25), annRow=me, annCol=cond, annColors=list("cond"=brewer.pal(3, "Set2")), hclustfun = function(d){hclust(d, method="ward.D2")}, distfun = "spearman",
    filename = sprintf("%s_DEHeatmap_FDR_%s.pdf", plot.name, p.thresh))

}
