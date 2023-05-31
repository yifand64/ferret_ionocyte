library(dplyr)
library(reshape)

#Based on Rebecca's code. 
# hypergeomtric test for changes in abundance between conditions.
do.celltype.abundance.plot <- function(celltype.id=NULL, condition.id=NULL, control.condition="Control", treatment.condition="Treatment", show.legend=T, ci=95,
                                      pval.cutoff.high=1e-5, pval.cutoff.low=0.05, ymax.expand=1,cols.use=NULL, nrow=NULL, celltype.order=NULL)
{
	if(is.null(celltype.id) | is.null(condition.id))
	{
		stop("Provide a celltype.id and condition.id")
	}  
	if(!(control.condition %in% condition.id) | !all(treatment.condition %in% condition.id))
	{
		stop("Provide a control and treatment conditions that are in the condition.id")
	}

	if(pval.cutoff.high > pval.cutoff.low){stop("pval.cutoff.high must be SMALLER than pval.cutoff.low")}
  
	data <- data.frame(condition = condition.id,cell.type = celltype.id,count=1)
	data.summarized <- data.frame(data %>% group_by(condition,cell.type) %>% summarise_each(funs(sum)))
	data.summarized <- cast(data.summarized,condition ~ cell.type)
	rownames(data.summarized) <- data.summarized$condition
	data.summarized <- data.summarized[,-1]
	data.summarized[is.na(data.summarized)] <- 0
	data.summarized$total <- rowSums(data.summarized) 
	data.summarized.fraction <- data.summarized/data.summarized$total

	info("Cell numbers in each group:")
	print(data.summarized)
	table.file = paste0(make.names(control.condition), "_vs_", make.names(paste(treatment.condition, collapse=", ")), "_numbers", ".txt")
	info(sprintf("Writing cell-abundance data to %s", table.file))
	write.table(data.summarized, file=table.file, sep="\t", quote=F)

	num.treatments = length(unique(treatment.condition))
	num.celltypes = length(unique(celltype.id))

	if(is.null(cols.use)){
		cols.use = colorRampPalette(brewer.pal(9, "Set1"))(num.celltypes)
	}
	print(sprintf("Measuring the effect of these %s treatments on the abundance of the %s celltypes", num.treatments, num.celltypes))

	enrichment.sig.pval <- matrix(NA, ncol=num.treatments,nrow=num.celltypes)
	depletion.sig.pval <- matrix(NA, ncol=num.treatments,nrow=num.celltypes)
	fold.change <- matrix(NA, ncol=num.treatments, nrow=num.celltypes)

	rownames(enrichment.sig.pval) <- unique(celltype.id)
	colnames(enrichment.sig.pval) <- treatment.condition
	rownames(depletion.sig.pval) <- unique(celltype.id)
	colnames(depletion.sig.pval) <- treatment.condition
	rownames(fold.change) <- unique(celltype.id)
	colnames(fold.change) <- treatment.condition

	for (curr.cell.type in unique(celltype.id)){
		for (exp in treatment.condition){
			# extract the numbers
		  	total.control.cells <- data.summarized[control.condition,'total']
			total.cells.curr.exp <- data.summarized[exp,'total']
			cells.of.curr.exp.in.curr.cell.type <- data.summarized[exp,curr.cell.type]
			cells.in.curr.cell.type <- sum(data.summarized[c(exp,control.condition),curr.cell.type])
			# test if any of the cells are missing from the current experiment or control
			if (data.summarized[control.condition,curr.cell.type] == 0 | cells.of.curr.exp.in.curr.cell.type == 0){
			  print(sprintf("This comparison of %s cells is omitted because the cell numbers for %s experiment are  %s and for %s are %s",
			                curr.cell.type, control.condition,data.summarized[control.condition,curr.cell.type],exp,cells.of.curr.exp.in.curr.cell.type))
			  next
			}
			# test for enrichment
			enrichment.sig.pval[curr.cell.type,exp] <- phyper(q=cells.of.curr.exp.in.curr.cell.type,k=cells.in.curr.cell.type,
			                  m=total.cells.curr.exp,n=total.control.cells,lower.tail=F)
			# test for depletion
			depletion.sig.pval[curr.cell.type,exp] <- phyper(q=cells.of.curr.exp.in.curr.cell.type,k=cells.in.curr.cell.type,
			                  m=total.cells.curr.exp,n=total.control.cells,lower.tail=T)
			# report fold change
			fold.change[curr.cell.type,exp] <- data.summarized.fraction[exp,curr.cell.type]/data.summarized.fraction[control.condition,curr.cell.type]
		}
	}

	data.for.barplot <- rbind(melt(-log10(enrichment.sig.pval)),melt(log10(depletion.sig.pval)))
	colnames(data.for.barplot) <- c('cell.type','experiment','pval')
	data.for.barplot$experiment <- factor(data.for.barplot$experiment, levels=treatment.condition)
	data.for.barplot$Color <- c(rep('enriched',num.celltypes*num.treatments),rep('depleted',num.celltypes*num.treatments))
	data.for.barplot <- na.omit(data.for.barplot)

	print(data.for.barplot)

	if(!is.null(celltype.order))
	{
		info("Reordering")
		data.for.barplot$cell.type = factor(data.for.barplot$cell.type, levels=celltype.order)
	}

	data.for.pval.tab = rbind(melt(enrichment.sig.pval), melt(depletion.sig.pval))
	print(head(data.for.pval.tab))
	colnames(data.for.pval.tab) <- c('cell.type','experiment','pval_raw')

	experiment = factor(data.for.pval.tab$experiment, levels=treatment.condition)
	print(table(experiment))
	data.for.pval.tab$experiment <- experiment
	data.for.pval.tab$direction <- c(rep('enriched',num.celltypes*num.treatments),rep('depleted',num.celltypes*num.treatments))
	data.for.pval.tab <- na.omit(data.for.pval.tab)
	table.file = paste0(make.names(control.condition), "_vs_", make.names(paste(treatment.condition, collapse=", ")), "_pvals", ".txt")
	info(sprintf("Writing p-value data to %s", table.file))
	write.table(data.for.pval.tab, file=table.file, sep="\t", quote=F, row.names=F)
	info("Drawing cp-value plot")
	P1 <- ggplot(data.for.barplot,aes(x=cell.type, y=pval)) + 
	  	geom_bar(aes(fill=Color),stat='identity') + 
	  	facet_wrap(~experiment, nrow=nrow) + 
	  	geom_hline(yintercept=c(-1.30103,1.30103)) +
	  	theme(panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border= element_blank(),
			panel.background =  element_rect(fill = 'white'),
			strip.background = element_blank(), 
			strip.text.x = element_text(size = 12, face="bold"),
			axis.text.x = element_text(angle=45, hjust=1),
			axis.title.x = element_text(vjust=-0.5),
			legend.key = element_blank(),
			text = element_text(size=12, colour="gray22"),
			axis.line.x = element_line(color="gray15", size = 0.75, lineend="round"),
			axis.line.y = element_line(color="gray15", size = 0.75, lineend="round")) +
	  	scale_fill_manual(values=c('enriched'='red','depleted'='navy'))

	info("Processing FC data")
	# make the same for fold change of relative abundance
	fold.change.ml <- melt(fold.change)
	colnames(fold.change.ml) <- c('cell.type','experiment','FC')
	fold.change.ml$experiment <- factor(fold.change.ml$experiment,levels=treatment.condition)
	fold.change.ml <- na.omit(fold.change.ml)

	info("Assigning significance stars")
	# make data frame for significance stars
	label.sig.twostar <- data.for.barplot[data.for.pval.tab$pval_raw < pval.cutoff.high,]
	if(nrow(label.sig.twostar)==0){
		warn("No very significant (**) changes!")
		there.are.twostar.changes = F
	}else{
		there.are.twostar.changes = T
		label.sig.twostar <- merge(label.sig.twostar,fold.change.ml,by=c('cell.type','experiment'))
		label.sig.twostar$FC <- label.sig.twostar$FC+0.05
		label.sig.twostar$label <- '**'
	}
	
	#print(label.sig)
	label.sig.onestar <- data.for.barplot[(data.for.pval.tab$pval_raw < pval.cutoff.low) & (data.for.pval.tab$pval_raw >= pval.cutoff.high),]
	if(nrow(label.sig.onestar)==0){
		warn("No significant (*) changes!")
		there.are.onestar.changes = F
	}else{
		label.sig.onestar$label <- '*'
		label.sig.onestar <- merge(label.sig.onestar,fold.change.ml,by=c('cell.type','experiment'))
		label.sig.onestar$FC <- label.sig.onestar$FC+0.05
		there.are.onestar.changes = T
	}
	

	if(!is.null(celltype.order))
	{
		info("Reordering")
		fold.change.ml$cell.type = factor(fold.change.ml$cell.type, levels=celltype.order)
	}

	print("Drawing FC plot")
	max_val = max(fold.change.ml$FC) + ymax.expand
	P2 <- ggplot(fold.change.ml,aes(x=cell.type,y=FC)) + 
	  	  geom_bar(aes(fill=cell.type),stat='identity') + 
	  	  facet_wrap(~experiment, nrow=nrow) + 
	  	  geom_hline(yintercept=c(1)) + xlab("") + 
	  	  # geom_hline(yintercept=c(2,0.5),linetype=2) +
	  	  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				panel.border= element_blank(),
				panel.background =  element_rect(fill = 'white'),
				strip.background = element_blank(), 
				strip.text.x = element_text(size = 12, face="bold"),
				axis.text.x = element_text(angle=45, hjust=1),
				axis.title.x = element_text(vjust=-0.5),
				legend.key = element_blank(),
				text = element_text(size=12, colour="gray22"),
				axis.line.x = element_line(color="gray15", size = 0.75, lineend="round"),
				axis.line.y = element_line(color="gray15", size = 0.75, lineend="round")) 
	

	if(there.are.twostar.changes){
		P2 = P2 + geom_text(data=label.sig.twostar,aes(label=label))
	}
	if(there.are.onestar.changes){
		P2 = P2 + geom_text(data=label.sig.onestar,aes(label=label))
	}

	if(show.legend){
		P2 <- P2 + scale_fill_manual('',values=cols.use, guide="legend")
	}else{
		P2 <- P2 + scale_fill_manual('',values=cols.use, guide="none")
	}
	
	return(list(P1,P2))
}


### test changes in celltype abundance using a Z-test relying
### on the central limit thereom. [ryan's idea]
test.abundance.ztest <- function(
	celltype.labels, 
	condition.labels, 
	batch.labels, 
	test.cell="Tuft", 
	compare=c("Control", "Day10 Helminth"), 
	plots=T,
	verbose=T)
{
	d = cbind.data.frame(celltype.labels, condition.labels, batch.labels)
	colnames(d) = c("celltype", "condition", "batch")
	if(!all(compare %in% d$condition)){stop(sprintf("'compare' arg must be levels of the condition labels! \n [%s]", paste(levels(d$condition), collapse=", ")))}
	d = d[d$condition %in% compare,]
	d$y = 0
	d$y[d$celltype==test.cell] <- 1
	d$batch = factor(d$batch)
	p = as.data.frame.matrix(table(d$y, d$batch)) 
	pc = as.data.frame.matrix(table(d$batch, d$condition)) 
	condition = colnames(pc)[unlist(apply(pc, 1, function(x){which(x>0)}))]
	
	#return(list("p"=p, "cond"=condition))
	print(p)
	print(pc)
	p = sweep(p, 2, colSums(p), "/")
	print(p)

	# for each celltype
	#look at notes of ryan's whiteboard
	#1. compute mouse-to-mouse variance (sigma_c)
	mu.ctl = unlist(p[2, condition==compare[1]])
	mu.trt = unlist(p[2, condition==compare[2]])
	
	print("mu.ctl:")
	print(mu.ctl)

	print("mu.trt:")
	print(mu.trt)

	var.ctl = var(mu.ctl)
	print("var.ctl")
	print(var.ctl)

	#2. compute s [s^2 = sigma_c^2/2]
	print("s")
	s = var.ctl / sqrt(2)
	print(s)

	#compute the Z-statistic
	z = (mean(mu.trt) - mean(mu.ctl)) / s
	p = 2*pnorm(-abs(z))

	info(sprintf("Z=%s, p=%s", z, p))
}


### conf.interval and point-estimate for proportion of a celltype, pooling
### estimates from multiple replicates (meta-analysis).
pool.ci = function(
	celltype.labels, 
	condition.labels, 
	batch.labels, 
	test.cell="Tuft"){}



test.abundance.poisson <- function(celltype.labels, 
	condition.labels, 
	batch.labels, 
	test.cell="Tuft", 
	verbose=F,
	coverage=NULL,
	coverage.agg=median,
	no.messages=F,
	neg.binom=F, #poisson if F, binom binomial if true
	compare=c("Control", "Day10 Helminth"))
{
	d = cbind.data.frame(celltype.labels, condition.labels, batch.labels)
	colnames(d) = c("celltype", "condition", "batch")
	if(!all(compare %in% d$condition)){stop(sprintf("'compare' arg must be levels of the condition labels! \n [%s]", paste(levels(d$condition), collapse=", ")))}
	if(!test.cell %in% d$celltype){stop(sprintf("'test.cell' arg must be a level of the celltype labels! \n [%s]", paste(levels(d$celltype), collapse=", ")))}
	d = d[d$condition %in% compare,]
	d$y = 0
	d$y[d$celltype==test.cell] <- 1
	d$batch = factor(d$batch)
	p = as.data.frame.matrix(table(d$y, d$batch)) 
	pc = as.data.frame.matrix(table(d$batch, d$condition)) 
	y = data.frame(t(p))
	colnames(y) = c("notCell", "isCell")
	y$total = rowSums(y[,1:2])
	y$condition = colnames(pc)[unlist(apply(pc, 1, function(x){which(x>0)}))]
	if(!is.null(coverage))
	{
		info(sprintf("Controlling for coverage [aggregate=%s]", as.character(substitute(coverage.agg))))
		cov = group.means(t(coverage), batch.labels)
		cov = data.frame(t(cov))
		colnames(cov) = "cov"
		cov = cov[rownames(y),]
		y$cov = cov
		if(neg.binom)
		{	
			library(MASS)
			psn = glm.nb(formula = isCell ~ condition + cov + offset(log(y$total)), data=y, maxit=1000)
		}else{
			psn = glm(formula = isCell ~ condition + cov, family = poisson, offset = log(y$total), data=y, maxit=1000)
		}
	}else{
		if(neg.binom)
		{	
			library(MASS)
			if(verbose){
				psn = glm.nb(formula = isCell ~ condition + offset(log(y$total)), data=y, maxit=1000, control=glm.control(trace = 3))# for lots of verbose output
			}else{
				psn = glm.nb(formula = isCell ~ condition + offset(log(y$total)), data=y, maxit=1000) 
			}
		}else{
			psn = glm(formula = isCell ~ condition, family = poisson, offset = log(y$total), data=y, maxit=1000)
		}
	}
	### make sure Poisson assumption is satisfied
	if(!neg.binom &! no.messages)
	{
		library(AER)
		d = dispersiontest(psn,trafo=1)
		if(d$p.value < 0.05) message(sprintf("%s cell numbers are likely overdispersed - p=%s, (Cameron & Trivedi, 1990) - better set neg.binom=T", test.cell, round(d$p.value,6)))
	}
	if(verbose)print(y)
	pv = coef(summary(psn))[,4]
	effect = coef(summary(psn))[,1]
	pv = pv[2]
	effect = effect[2]
	return(list("poisson"=psn, "p.value"=pv, 
		"effect.estimate"=exp(effect), 
		"conf.95"=exp(confint.default(psn)), 
		"conf.99"=exp(confint.default(psn,level = 0.99))))
}

### same function. coded in a more clear way. (not sure what was going on the first time)
test.abundance.logit <- function(seur, 
	celltype.label, 
	condition.label,
	batch.label=NULL,
	use.nGene=F,
	use.nUMI=F,
	test.cell="Tuft", 

	compare=c("Control", "Day10 Helminth"))
{
	library(lme4)
	seur@data.info$nUMI = colSums(seur@raw.data)
	info("Fetching data from seurat object")
	m = fetch.data(seur, c(condition.label, batch.label, celltype.label, "nGene", "nUMI"))
	m = m[m[, condition.label] %in% compare,]
	colnames(m) = c("condition", "batch", "celltype", "nGene", "nUMI")
	m$nGene = m$nGene / max(m$nGene) # use CDR like MAST
	m$y = 0
	m$y[m$celltype==test.cell] <- 1
	mc = m[m$Condition %in% compare,]
	
	info(sprintf("Fitting model [use nGene=%s]", use.nGene))
	if(use.nGene)
	{
		##glmer parameters # http://stats.idre.ucla.edu/r/dae/mixed-effects-logistic-regression/
		logit = glmer(y ~ condition + nGene + (1|batch), data = m, 
			family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
	}else{
		##glmer parameters # http://stats.idre.ucla.edu/r/dae/mixed-effects-logistic-regression/
		logit = glmer(y ~ condition + (1|batch), data = m, family = binomial, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
	}
	pv = coef(summary(logit))[,4]
	effect = coef(summary(logit))[,1]
	pv = pv[2]
	effect = effect[2]
	return(list("logit"=logit, "p.value"=pv, "effect.estimate"=effect))
}





# ### test changes in celltype abundance using a logistic regression to
# ### control for mouse-to-mouse variation (batch). we convert the cell labels
# ### to a table of cell-level observations, in which each observed cell either is (Y=1)
# ### or is not (Y=0) the celltype of interest. We can then use a logistic
# ### regression to predict Y
# test.abundance.logit <- function(
# 	celltype.labels, 
# 	condition.labels, 
# 	batch.labels, 
# 	complexity=NULL,
# 	test.cell="Tuft", 
# 	compare=c("Control", "Day10 Helminth"), 
# 	plots=T,
# 	verbose=T)
# {
# 	m =  as.data.frame.matrix(table(batch.labels, celltype.labels))
# 	m$other = rowSums(m[, colnames(m) != test.cell])
# 	p = as.data.frame.matrix(table(batch.labels, condition.labels)) 
# 	m$condition = colnames(p)[unlist(apply(p, 1, function(x){which(x>0)}))]
# 	m = m[m$condition %in% compare, c(test.cell, "other", "condition")]
# 	m$batch = rownames(m)
# 	if(verbose){
# 		print(m)
# 	}
# 	rownames(m) <- NULL
# 	m = melt(m)
# 	m <- m[rep(sequence(nrow(m)), m[["value"]]), ]
# 	m[["value"]] <- NULL
# 	m$y = mapvalues(m$variable, c(test.cell, "other"), c(1, 0))
# 	print(head(m))
# 	if(!is.null(complexity))
# 	{
# 		m$complexity = complexity
# 	}
	
# 	library(lme4)
# 	### 
# 	## fit a variable intercept multilevel logistic regression.
# 	## the 'random effects' term (1|batch) tells the solver to
# 	## fit a variable intercept for each mouse (batch).
	
	
	
# 	if(!is.null(complexity))
# 	{
# 		logit = glmer(y ~ condition + complexity + (1|batch), data = m, family = binomial(link='logit'))
# 	}else{
# 		logit = glmer(y ~ condition + (1|batch), data = m, family = binomial(link='logit'))
# 	}
	
# 	if(plots){
# 		d = coef(logit)$batch
# 		d = cbind(d, subset(as.data.frame(table(m$batch)), Var1 %in% rownames(d)))
# 		colnames(d) = c("intercept", "slope", "mouse", "observations")
# 		p1 = ggplot(d, aes(x=observations, y=intercept, label=mouse)) + geom_point() + geom_text(vjust=2)
# 		ggsave(p1, filename="variable_intercepts.pdf")
# 	}
	

# 	if(verbose){
# 		print(summary(logit))
# 		#print(confint(logit)) ## takes a while
# 	}
# 	#alias(logit)
# 	pv = coef(summary(logit))[,4]
# 	effect = coef(summary(logit))[,1]
# 	pv = pv[2]
# 	effect = effect[2]
# 	return(list("logit"=logit, "p.value"=pv, "effect.estimate"=effect))
# }


### Set up 2x2xk contingency table
### example use:
### m = conting(all@data.info$Graph_K500_named, all@data.info$Condition, all@data.info$Batch, test.cell = "Goblet", compare = c("Control", "Day3 Helminth"))


### literature:
# https://en.wikipedia.org/wiki/Cochran%E2%80%93Mantel%E2%80%93Haenszel_statistics#cite_note-agresti-1
# http://mathdept.iut.ac.ir/sites/mathdept.iut.ac.ir/files/AGRESTI.PDF [p231--237]
# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/mantelhaen.test.html

### consider trying to simplify the construction of the contingency table using 'xtabs' 
# http://stat.ethz.ch/R-manual/R-devel/library/stats/html/xtabs.html
test.abundance.mantelhaen <- function(
	celltype.labels, 
	condition.labels, 
	batch.labels, 
	test.cell="Tuft", 
	compare=c("Control", "Day10 Helminth"), 
	verbose=T)
{
	library(abind)
	m =  as.data.frame.matrix(table(batch.labels, celltype.labels))
	m$other = rowSums(m[, colnames(m) != test.cell])
	p = as.data.frame.matrix(table(batch.labels, condition.labels)) 
	m$condition = colnames(p)[unlist(apply(p, 1, function(x){which(x>0)}))]
	m = m[m$condition %in% compare, c(test.cell, "other", "condition")]
	#print(m)

	keep.pair.if.right.comparison = function(pair){
		tab = m[pair,]
		if(length(unique(tab$condition)) == 1){
			return(NA)
		}else{
			cond.col = unlist(tab$condition)
			# print(cond.col)
			num.in.order = sum(cond.col == compare)
			# print(num.in.order)
			if(num.in.order==2){
				return(tab[,1:2])
			}else{
				return(tab[c(2,1),1:2])
			}
		}
	}
	tabs = apply(data.frame(combn(1:nrow(m),2)), 2, keep.pair.if.right.comparison)
	tabs = tabs[!is.na(tabs)]
	comps = unlist(lapply(tabs, function(x){paste(rownames(x), collapse="_")}))
	#print(comps)
	#dimnames = list(condition=compare, celltype=c(test.cell, "other"), comparison=comps)
	names(tabs) = comps
	arr = do.call(abind, c(tabs, along = 3))

	if(verbose){print(arr)}

	if(verbose){info("Conditional odds ratios [assumed equal]:")}
	woolf.pval = woolf(arr)
	if(verbose){
		if(woolf.pval > 0.05){info(sprintf("Woolf test for interactions: %s", woolf.pval))}else{
			warn(sprintf("Woolf test for interactions: %s! \nassumption of homogeneous (conditional) odds ratios may be violated!", woolf.pval))
		}
		print(apply(arr, 3, function(x) (x[1,1]*x[2,2])/(x[1,2]*x[2,1])))
	}
	mantelhaen.test(arr)
}



#the Woolf test for interaction:
woolf <- function(x) {
  x <- x + 1 / 2
  k <- dim(x)[3]
  or <- apply(x, 3, function(x) (x[1,1]*x[2,2])/(x[1,2]*x[2,1]))
  w <-  apply(x, 3, function(x) 1 / sum(1 / x))
  1 - pchisq(sum(w * (log(or) - weighted.mean(log(or), w)) ^ 2), k - 1)
}



