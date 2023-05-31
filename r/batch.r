

# calculate the batch composition of each cluster, and draw a set of pie charts
batch.composition <- function(clusters, batches, 
              title="Batch composition of clusters",
              get.freq=T,
              ret.data=T,
              pie=F,
              xlab="",
              nrow=NULL,
              cex.title=10,
              cex.label=8, 
              legend.title="Batch",
              pdf=NULL, 
              label.frac=T,
              sort.alphabetical=T,
              sort.by = NULL,
              cols.use=NULL,
              order.groups=NULL,
              width=8, 
              height=7)
{
    df = NULL
    for(cl in unique(clusters))
    {
      cells = which(clusters==cl)
      b = batches[cells]
      n = length(cells)
      if(get.freq){
        d = data.frame(table(b)/ n)
      }else{
        d = data.frame(table(b))
      } 
      d$cn = rep(cl, nrow(d))
      d$cluster.size = rep(n, nrow(d))
      d$cluster.name = rep(paste(cl, sprintf("\n[n=%s]", n)), nrow(d))
      if(is.null(df))
      {
        df = d
      }else
      {
        df = rbind(df, d)
      }
    }
    df$Frac = df$Freq/df$cluster.size
    # put the batches in alphabetical order
    if(sort.alphabetical){
        new.levels = sort(unique(as.character(unlist(df$b))))
    }else{
        new.levels = levels(df$b)
    }
    if(!is.null(order.groups))
    {
        if(all(order.groups %in% unique(df$cn)))
        {
            info("Re-ordering groups")
            # grab the correct new group labels, including the [n=whatever] part
            use.order.groups = unique(df$cluster.name)[match(order.groups, unique(df$cn))]
            df$cluster.name = factor(as.character(unlist(df$cluster.name)), levels=use.order.groups)
        }else{
              warn(sprintf("Could not re-order groups. %s are not part of the groups", 
              paste( order.groups[!order.groups %in% unique(df$cluster.name)], collapse=", ")))
        }
    }
    df$Batch = factor(as.character(unlist(df$b)), levels=new.levels)
    n.cols = length(new.levels)
    ### calculate the position of text labels:'
    library(dplyr)
    df <- ddply(df, .(cn), 
       transform, pos = cumsum(Freq) - (0.5 * Freq)
    )
    if(ret.data){
      return(df)
    }
    if(pie){
        p <- ggplot(df, aes(x=1, y=Freq, fill=Batch)) +
          ggtitle(title) +
          theme_bw() + 
          theme(strip.text.x = element_text(size = cex.title, face="bold"), 
              panel.border= element_blank(),
              panel.background =  element_blank(),
              strip.background = element_blank()) + 
          scale_fill_manual(legend.title, values=colorRampPalette(brewer.pal(n.cols, "Set1"))(n.cols)) + 
          coord_polar(theta='y') + 
          theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x=element_blank(),
            panel.grid  = element_blank()) + 
              # black border around pie slices
             geom_bar(stat="identity", color='black') +
            # remove black diagonal line from legend
            guides(fill=guide_legend(override.aes=list(colour=NA))) + facet_wrap(~cluster.name, nrow=nrow) 
    }else{
       print(legend.title)
        p <- ggplot(df, aes(x=1, y=Freq, fill=Batch, label=Freq)) +
          ggtitle(title) +
          theme_bw() + 
          theme(strip.text.x = element_text(size = cex.title, face="bold"),
              panel.border= element_blank(),
              panel.background =  element_blank(),
              strip.background = element_blank() ) + 
          scale_fill_manual(legend.title, values=colorRampPalette(brewer.pal(n.cols, "Set1"))(n.cols)) + 
          theme(axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.title.x=element_blank(),
                legend.key = element_blank(),
            # axis.line.x = element_line(color="gray15", size = 0.75, lineend="round"),
            # axis.line.y = element_line(color="gray15", size = 0.75, lineend="round"),
            panel.grid  = element_blank()) + 
            # black border around bars, but not in the legend.
            geom_bar(stat="identity", position="fill") +
            geom_bar(stat="identity", color='black', show.legend=FALSE, position="fill") + facet_wrap(~cluster.name, nrow=nrow) 

          if(label.frac)
          {
              library(scales)
              p = p +  geom_text(aes(label = percent(signif(Freq,3))), size=cex.label, position = position_fill(vjust = 0.5))
          }
          p + ylab(xlab)
    }
      if(!is.null(cols.use))
      {
        # print(cols.use)
        # print(order(table(df$Batch)))
        # cols.use = cols.use[order(table(df$Batch))]
        # print(cols.use)
        p = p + scale_fill_manual(legend.title, values=cols.use) 
      }
      p = p + xlab(xlab)
    if(!is.null(pdf))
    {
      ggsave(p, filename=pdf, width=width, height=height)
    }else
    {
      print(p)
    }
}


scatter.plot <- function(data, x.col, y.col, show.points=T, show.density=T,colour.pts =F, use.smooth=F, pt.size=1,
                          pdf=NULL, width=9, height=8, cols=rich12equal, xlab=NULL, ylab=NULL, show.outliers=T, use.image=T,
                          alpha = 1, interpolate=T, xlim=c(0,10), ylim=c(0,10),
                          upper = 1, lower = -1, use.alpha=F)
{
    if(is.null(xlab)){xlab=x.col}
    if(is.null(ylab)){ylab=y.col}

    x = data[, x.col]
    y = data[, y.col]
    if(use.smooth)
    {
        info(sprintf("use.smooth is: %s", use.smooth))
        if(!is.null(pdf)){pdf(pdf, width=width,height=height)}
          info(sprintf("saving to: %s", pdf))
          smoothScatter(x,y, nrpoints=Inf, cex=pt.size, colramp=colorRampPalette(cols), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim) + 
            text(2, 7, cex=2, sprintf("r = %s", signif(cor(x,y), 3)))
        if(!is.null(pdf)){dev.off()}
        return(NULL)
    }
    data$label = rownames(data)
    data[, "density"] <- densCols(x,y, colramp=colorRampPalette(cols))

    print(head(data))
    g = ggplot(data, aes_string(x=x.col, y=y.col, label="label"))
    
    if(show.points)
    {
        info("Drawing points")
        # draw a scatter plot to visualise the correlation. colour pts by density
         if(colour.pts){
             info("Colouring points by density")
             g = g + geom_point(aes(color=density), size=pt.size) + 
              scale_color_identity(guide="colorbar")  
         }else{
             g = g + geom_point(size=pt.size) 
          }
        
    }
    if(show.density){
        info("Drawing density")
        library(MASS)
        dens <- kde2d(x, y, n=200)
        dfdens <- data.frame(expand.grid(x=dens$x, y=dens$y), z=as.vector(dens$z))

        if(use.image){
            #image(dens, col=dens.cols)
            info(sprintf("Alpha=%s", alpha))
            g = g + geom_raster(data=dfdens, aes(x=x, y=y, fill=z), alpha=alpha, interpolate=interpolate) 
        }else{
            g = g+ stat_contour(aes(x=x, y=y, z=z, fill=..level.., alpha=..level..), 
            data= dfdens, geom="polygon")
        }
        
        g = g + scale_fill_gradientn("Density", colors=cols, guide="colorbar", labels=NULL) +
            guides(alpha=FALSE) + theme_bw() 
        return(g)
        # otherwise density near 0,0 overwhelms the plot
        #g = g + scale_x_continuous(limits=c(1,11)) +  scale_y_continuous(limits=c(1,11)) 
    }

    

    if(show.outliers)
    {
        
        low.ol = data[y < x + lower,]
        high.ol = data[y > x + upper,]

        info(sprintf("Showing %s high and %s low outlier genes", nrow(high.ol), nrow(low.ol)))

        g = g + geom_abline(intercept = upper, slope = 1, linetype="dotted") + 
          geom_abline(intercept = lower, slope = 1, linetype="dotted") + theme_bw() + 
          
          #label points
          geom_text(data=high.ol, size=2, vjust=1.5) + 
          geom_text(data=low.ol, size=2, vjust=1.5) + 
          # colour points
          geom_point(data=high.ol, colour="red3") + 
          geom_point(data=low.ol, colour="dodgerblue4") 
          
    }
    g = g + theme_bw() +  
          annotate("text", x=6, y=1, label= sprintf("R=%s", signif(cor(x,y), 3)), fontface="bold") + 
          xlab(xlab) +  ylab(ylab) + theme(axis.line.x = element_line(color="gray15", size = 0.75, lineend="round"),
                                              axis.line.y = element_line(color="gray15", size = 0.75, lineend="round")) + 
          theme(
                panel.border= element_blank(),
                axis.title.x = element_text(vjust=-0.5),
                #plot.margin = unit(c(1.2,1.2,1.2,1.2), "cm"),
                text = element_text(size=20, colour="gray22"))
    if(!is.null(pdf))
    {
        ggsave(g, filename=pdf, width=width, height=height)
    }
    return(g)
}



# look at the correlation of average expression between
# batches of a single-cell (seurat) object
batch.scatter <- function(s, data.is.tpm=T, show.pop.miseq=F, expression.bins=c(0, 1, 2, 4, 8))
{
    library(Hmisc)
    library(GGally)

    if(data.is.tpm){
        data = log2(1+s@raw.data)
    }else{
        data = log2(1+tpm(s@raw.data))   
    }
    gm = data.frame(group.means(data, groups = s@data.info$Batch))

    if(show.pop.miseq){
        pop.miseq = read.delim("/Users/ahaber/Desktop/projects/GutCircuits/ahaber/FROM_USERS_AHABER/data/population_cell_types_miseq/data/pipeline_output/miseq/samples/merged_counts_table.txt")
        gm$GENE_SYMBOL = rownames(gm)
        gmm = merge(gm, pop.miseq, by = "GENE_SYMBOL")
        gmm$LogTPM_Bulk_Miseq = log2(1+gmm$TPM_bulk)
        gm = gmm[, c(1:11, ncol(gmm))]
        gm$Expr.Bin = cut2(g$LogTPM_Bulk_Miseq, cuts=expression.bins)
    }else{
        gm$Expr.Bin = cut2(rowMeans(data), cuts=expression.bins)
    }
    
    ggpairs(gm, showStrips = T, mapping = aes(colour=Expr.Bin), 
            lower = list(continuous = wrap("points", alpha = 0.3,    size=0.1), 
                         combo = "box"))
}


# ### the commented code is my original attempt.
# ### the rest is based on code Karthik sent over.
batch.normalise.comBat <- function(counts, batch.groups, max.val=0, par.prior=T, preserve.zeros=T) #,filter.low.variance.genes=F)
{
    
    ##### [Adam added this to deal with zero variance rows]
    info("Computing variances")
    v = RowVar(counts)
    if(any(v == 0)){
        zero_var_genes = names(v)[v == 0]
        warn(sprintf("%s genes had zero variance and could not be corrected. left alone", length(zero_var_genes)))
        print(zero_var_genes)
        uncorrected_rows = T
    }else{
        uncorrected_rows = F
    }
    #####

    # if(filter.low.variance.genes)
    # {
    #     n = nrow(counts)
    #     info("Filtering low-variance genes is ON")
    #     batch.vars = apply(t(counts),2,function(x) tapply(x,batch.groups,var)) 
    #     batch.vars = t(na.omit(batch.vars))
    #     print(head(batch.vars))
    #     use.genes = which(rowSums(batch.vars>0)==ncol(batch.vars)) ##non-zero variance in all batches
    #     counts = counts[use.genes,]
    #     info(sprintf("Removed %s low-variance genes, %s remain", n-nrow(counts), nrow(counts)))
    # }

    library(sva)
    batch.groups = factor(batch.groups) ## drop zero levels
    batch.id = 1:length(unique(batch.groups))
    names(batch.id) = unique(batch.groups)
    batch.ids = batch.id[batch.groups]
    
    correct.data = data.frame(matrix(0L, nrow = dim(counts)[1], ncol = dim(counts)[2]))      
    correct.data[v > 0, ] = ComBat(counts[v > 0, ],batch.ids, prior.plots=FALSE, par.prior=par.prior)
    if(uncorrected_rows){correct.data[v == 0,] <- counts[v == 0,]}
    if(max.val > 0){correct.data[correct.data > max.val] = max.val}
    if(preserve.zeros){correct.data[counts==0] <- 0}
    return(as.data.frame(correct.data))
    
}


RowVar <- function(x, ...) {
    library(Matrix)
    Matrix::rowSums((x - Matrix::rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

within.celltype.batchcorrect <- function(counts, celltype.labels, batch.labels, 
  use.combat=T,
  max.val=6, 
  par.prior=T)
{
    
    info(sprintf("Use comBat -->  %s", use.combat))
    groups = as.character(unique(celltype.labels))
    rval = data.frame(matrix(0L, nrow = dim(counts)[1], ncol = dim(counts)[2]))
    for(g in groups)
    {
        ingroup = which(celltype.labels==g)
        raw = counts[,ingroup]
        bl = batch.labels[ingroup]
        info(sprintf("Correcting %s %s cells", ncol(raw), g))
        info("Batch labels: ")
        print(table(bl))
        if(use.combat){
          corrected = batch.normalise.comBat(raw, bl, par.prior=par.prior)
        }else{
          corrected = batch.normalise(raw, bl, par.prior=par.prior)
        }
        
        rval[,ingroup] <- corrected
        cat("\n")
    }
    rval
}

# follow Carl De Boer's batch normalisation procedure.
# This works slightly differently in log space, where it is an additive scale 
# rather than a multiplicative scale - basically, I shift the centers of each batch 
# towards each other for each gene by subtracting half the difference of means 
# (in log(TPM/10+1) space) from one batch and adding it to the other.  This has 
# the effect of making the distributions of gene expression for each gene have the same mean.  
# This is similar to just removing the PC, but is better in my opinion because we remove the 
# known confounder directly rather than having PCA find it for us (which might include 
# other things as well and might not perfectly correspond to the batches).
batch.normalise <- function(counts, batch.groups, ignore.small=10, fn=mean, preserve.zeros=F)  
{
  info("Running batch normalisation")
  print(table(batch.groups))
  info("Finding batch means")
  batch.averages <- aggregate(t(counts), by=list(batch.groups), FUN=fn, na.rm=TRUE)
  rownames(batch.averages) = batch.averages$Group.1
  batch.averages$Group.1 = NULL
  batch.averages = data.frame(t(batch.averages), check.names=F)
  info("Top few batch means:")
  print(head(batch.averages))

  info("Finding correction")
  
  if(ignore.small > 0){

      large.enough = which(table(batch.groups) >= ignore.small)
      too.small = which(table(batch.groups) < ignore.small)
      means = rowMeans(batch.averages[, large.enough])
      if(length(large.enough) < length(unique(batch.groups)))
      {
          warn(sprintf("Following groups were not used in the correction\nas they had less than %s samples", ignore.small))
          print(table(batch.groups)[too.small])
      }
  }else{
      means = rowMeans(batch.averages)
  }
  correction = batch.averages - means
  cor.order = order(correction, decreasing=T)
  info("Largest correction")
  print(head(batch.averages[cor.order,]))

  adjusted.averages = batch.averages - correction
  info("Applying correction")
  corrected.cts = counts
  for(group in colnames(correction))
  {
      corrected.cts[, batch.groups==group] = counts[, batch.groups==group] - correction[, group]
  }
  
  if(preserve.zeros){corrected.cts[counts==0] <- 0}
  return(corrected.cts)
}




