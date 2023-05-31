library(igraph)
library(edgebundleR)
library(ggraph)


edge_plot = function( genes, 
                      celltype = NULL, 
                      ix_path="~/ref/ix/combined_interactions.txt", 
                      organism = "mouse", 
                      width=600,
                      fontsize=10,
                      use.cols=NULL){
  
  if(!file.exists(ix_path)){
    stop('ref folder does not exist in home directory')
  }
  interactions = read.table(ix_path, header = T, stringsAsFactors = FALSE)
  if(organism == 'mouse'){
    interactions$gene1 = interactions$gene1_mouse
    interactions$gene2 = interactions$gene2_mouse
  }
  else if(organism == 'human'){
    interactions$gene1 = interactions$gene1_human
    interactions$gene2 = interactions$gene2_human
  }
  else{
    stop("Unrecognized organim argument. Use mouse or human.")
  }
  dex = interactions$gene1 %in% genes + interactions$gene2 %in% genes
  double_interactions = interactions[which(dex ==2), ]
  double_interactions$gene1 = as.character(double_interactions$gene1)
  double_interactions$gene2 = as.character(double_interactions$gene2)
  rel2 = double_interactions[,c("gene1", "gene2")]
  if(!is.null(celltype)){
    d = as.data.frame(cbind(as.character(genes), as.character(celltype)))
    colnames(d) = c("gene", "cell")
    rownames(d) = as.character(d$gene)
    double_interactions$cell_type_1 = d[as.character(double_interactions$gene1), "cell"]
    double_interactions$cell_type_2 = d[as.character(double_interactions$gene2), "cell"]
    m = d %>% filter(gene %in% c(double_interactions$gene1, double_interactions$gene2))
    m$cell =gsub("\\.", "-", m$cell)
    if(is.null(cell_order)){cell_order=as.character(unique(m$cell))}else{cell_order = gsub("\\.", "-", cell_order)}
    d2 <- structure(list(ID = paste0(m$cell, ".", m$gene), Loc =factor(m$cell, levels = cell_order)), .Names = c("ID", "Loc"), class = "data.frame",row.names = c(NA, -length(m$gene)))
    d2$key = factor(m$gene, levels = m$gene)
    rel2$gene1 = paste0(gsub("\\.", "-", double_interactions$cell_type_1),".",double_interactions$gene1)
    rel2$gene2 = paste0(gsub("\\.", "-", double_interactions$cell_type_2),".",double_interactions$gene2)
    g2 <- graph.data.frame(rel2, directed=F, vertices=d2)
    clr <- factor(V(g2)$Loc)
    if(is.null(use.cols)){use.cols=distinct.cols(n=length(unique(m$cell)))}
    levels(clr) <- use.cols
    V(g2)$color <- as.character(clr)
    edgebundle( g2, width=width, fontsize = fontsize)->eb
    return(eb)
  }else{
    g2 <- graph.data.frame(rel2, directed=F)
    edgebundle( g2, width=width, fontsize = fontsize)->eb
    return(eb)
  }
}

node_plot = function(genes, celltype, cell_nums, ix_path="~/ref/ix/combined_interactions.txt", organism = "mouse", use.cols=NULL, node_size = c(5,25), text_size=8,
                     edge_color="gray50", edge_width = c(1,4), label_line_color = "transparent"){
  
  if(!file.exists(ix_path)){
    stop('ref folder does not exist in home directory')
  }
  
  interactions = read.table(ix_path, header = T, stringsAsFactors = FALSE)
  
  if(organism == 'mouse'){
    interactions$gene1 = interactions$gene1_mouse
    interactions$gene2 = interactions$gene2_mouse
  }
  
  else if(organism == 'human'){
    interactions$gene1 = interactions$gene1_human
    interactions$gene2 = interactions$gene2_human
  }
  
  else{
    stop("Unrecognized organim argument. Use mouse or human.")
  }
  
  dex = interactions$gene1 %in% genes + interactions$gene2 %in% genes
  double_interactions = interactions[which(dex ==2), ]
  double_interactions$gene1 = as.character(double_interactions$gene1)
  double_interactions$gene2 = as.character(double_interactions$gene2)
  
  d = as.data.frame(cbind(as.character(genes), as.character(celltype)))
  colnames(d) = c("gene", "cell")
  rownames(d) = as.character(d$gene)
  
  double_interactions$cell_type_1 = gsub("\\.", "-", d[as.character(double_interactions$gene1), "cell"])
  double_interactions$cell_type_2 = gsub("\\.", "-", d[as.character(double_interactions$gene2), "cell"])

  
  g <- graph(as.vector(t(double_interactions[,c("cell_type_1", "cell_type_2")])))
  edge_attr(g, "label") = as.character(double_interactions$interacting_pair)
  
  g2 = simplify(g)
  
  E(g2)$weight = sapply(E(g2), function(e){length(all_shortest_paths(g, from=ends(g2, e)[1], to=ends(g2, e)[2])$res)})
  
  V(g2)$n_cells = cell_nums[V(g2)$name]
  
  if(is.null(use.cols)){use.cols=distinct.cols(n=length(V(g2)))}
  
  gp = ggraph(g2, layout = "fr") + 
    geom_edge_link(alpha=.9, color=edge_color, aes(width=weight)) + 
    geom_node_point(aes(size=n_cells), color=use.cols) + scale_size(range = node_size) + 
    geom_node_text(aes(label=name), size=text_size, repel = T, force=20, segment.color=label_line_color) +
    theme(legend.position = "bottom") + 
    scale_edge_width(range = edge_width) 
  
  return(gp)
  
}
