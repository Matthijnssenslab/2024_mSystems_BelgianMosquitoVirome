#' plot_relabund
plot_relabund <- function (physeq, x = "Sample", y = "Abundance", fill = NULL,
                          color=NULL,
                          title = NULL, facet_grid = NULL, facet_wrap = NULL) 
{
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill, color=color))
  p = p + geom_bar(stat = "identity", position = "stack")
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(facet_wrap)) {
    p <- p + facet_wrap(facet_wrap, scales = "free_x")
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

#' heatmap function
plot_abundance <- function (obj, n, norm = TRUE, log = TRUE, fun = sd, rowclust = "bray", colclust = "bray", ...) 
{
  mat <- returnAppropriateObj(obj, norm, log)
  otusToKeep <- which(rowSums(mat) > 0)
  otuStats <- apply(mat[otusToKeep, ], 1, fun)
  otuIndices <- otusToKeep[order(otuStats, decreasing = TRUE)[1:n]]
  mat2 <- mat[otuIndices, ]
  mat3 <- mat2[order(row.names(mat2)),]
  mat4 <- as.data.frame(t(mat3))
  Heatmap(mat3, clustering_distance_columns = vegdist(mat4, method = colclust), clustering_distance_rows = vegdist(mat3, method = rowclust), ...)
}

#' Change labels of scientific format
scientific_10 <- function(x) {
  parse(text=gsub("1e\\+", "10^", scales::scientific_format()(x)))
}

#' helper function str_split
subset_safely <- function(x, index) {
  if (length(x) < index) {
    return(NA_character_)
  }
  x[[index]]
}

#' str_split
str_split_n <- function(string, pattern, n) {
  out <- str_split(string, pattern)
  vapply(out, subset_safely, character(1L), index = n)
}

#' function to clean the metadata of phylogenetic trees
clean_phylometa <- function(tree, metadata){
  namesvector <- c()
  for (i in 1:length(tree@phylo$tip.label)){
   #print(str_split_n(tree@phylo$tip.label[i], "\\|", 1))
   strname <- str_split_n(tree@phylo$tip.label[i], "\\|", 1)
   namesvector[i] <- strname
  }
  metadata <- metadata[is.element(metadata$Accession, namesvector),] #changed %in%

  metadata <- metadata %>% 
    mutate(tipname=paste(Accession, Species, sep="|")) %>% 
    #replace_na(list(Genus="unclassified_g", Family="unclassified_f")) %>%
    replace_na(list(Genus="unclassified", Family="unclassified")) %>%
    replace(is.na(.), "unknown") %>% 
    mutate_all(list(~str_replace_all(., ":.*", "")))
  
  return(metadata)
}

clean_phylometa_negev <- function(tree, metadata){
  namesvector <- c()
  for (i in 1:length(tree@phylo$tip.label)){
    #print(str_split_n(tree@phylo$tip.label[i], "\\|", 1))
    strname <- str_split_n(tree@phylo$tip.label[i], "_\\|", 1)
    namesvector[i] <- strname
  }
  metadata <- metadata[is.element(metadata$Accession, namesvector),] #changed %in%
  
  metadata <- metadata %>% 
    mutate(tipname=paste(Accession, Species, sep="|")) %>% 
    #replace_na(list(Genus="unclassified_g", Family="unclassified_f")) %>%
    replace_na(list(Genus="unclassified", Family="unclassified")) %>%
    replace(is.na(.), "unknown") %>% 
    mutate_all(list(~str_replace_all(., ":.*", "")))
  
  return(metadata)
}

#' Function to group the data in a tree based on value in the metadata 
group_phylotaxa <- function(tree, metadata, group="Family"){
  
  tree@phylo$tip.label <- gsub(tree@phylo$tip.label, pattern = "_", replacement=" ")
  tree@phylo$tip.label <- gsub(tree@phylo$tip.label, pattern = " \\|", replacement="\\|")
  tree@phylo$tip.label <- gsub(tree@phylo$tip.label, pattern = "P ", replacement="P_", perl=T)
  tree@phylo$tip.label <- gsub(tree@phylo$tip.label, pattern = "P_virus", replacement="P virus", perl=T)
  tree@phylo$tip.label <- gsub(tree@phylo$tip.label, pattern = "Parry s", replacement="Parry's", perl = T)

  for (g in group){
    cls=list()
    cls=list("NODE"=tree@phylo$tip.label[grepl("NODE", tree@phylo$tip.label, fixed = T)])
    for (i in metadata[[g]]){
      cls[i] <- list(metadata[metadata[[g]]==i,]$tipname) 
    } 
    sublist <- cls[names(cls) %in% c(unique(metadata[[g]]), "NODE")]
    tree <- groupOTU(tree, sublist, group=g)
  }
  return(tree)
}

#' Function to get the MRCA node in a tree of specified taxonomy level
#' Best used with ICTV confirmed species
mrca_ictv <- function(tree, meta, group, subset = NULL){
  treeplot <- ggtree(tree)
  tipdata <- subset(treeplot$data, is.element(treeplot$data$label, meta$tipname))

  if(!is.null(subset)){
    #ugroup <- unique(tipdata[[group]])[!is.element(unique(tipdata[[group]]), subset)]
    ugroup <- names(table(tipdata[[group]])[table(tipdata[[group]])>1])
    ugroup <- ugroup[!is.element(ugroup, subset)]
  } else {
    #ugroup <- unique(tipdata[[group]])
    ugroup <- names(table(tipdata[[group]])[table(tipdata[[group]])>1])
  }
  
  mat <- matrix(nrow = length(ugroup), ncol=2)
  for (i in 1:length(ugroup)){
    mat[i,] <- c(getMRCA(tree@phylo, tipdata[tipdata[[group]] == ugroup[i],]$node), as.character(ugroup[i]))
  }
  df <- as.data.frame(mat)
  colnames(df) <- c("node", "id")
  #df <- df[!is.na(as.numeric(as.character(df$node))),]
  rownames(df) <- NULL
  df <- transform(df, node = as.numeric(node))
  return(df)
}

#' Plot phylogenetic trees functions
plot_phylotree <- function(tree, col, shape, cladedf=NULL, labels, plainlabels=NULL, branchwidth=.25, 
                           anno.size=1.5, ladderize=T, support=70, scalewidth=0.5,
                           label=label, group="Family", align=T)
  {
  p <- ggtree(tree, ladderize = ladderize, size=branchwidth)+
    geom_tippoint(aes(color=!!as.name(group), shape=!!as.name("Classification")), size=.7)+
    geom_nodelab(aes(label=as.numeric(label), subset = !is.na(as.numeric(label)) & as.numeric(label) < 90), 
                 size=anno.size, hjust=-.01) +
    #geom_nodepoint(aes(fill=as.numeric(label), subset = !is.na(as.numeric(label))), 
    #               shape=23, color="transparent", size=anno.size-.5)+
    scale_color_manual(values = col, 
                       #name="NCBI classification",
                       name="",
                       breaks = labels,
                       labels=toexpr(c(labels[labels!="NODE"],"Belgian mosquitoes"), 
                                     plain=plainlabels))+
    scale_shape_manual(values = shape, 
                       name="")
  if (!is.null(cladedf)) {
    p <- p + 
      geom_cladelab(data=cladedf, mapping=aes(node=node, label=id), color="black",
                           fontface="italic", fontsize=2, barsize=.25, align=align, offset=.1)#+
      #geom_nodepoint(aes(subset= node %in% cladedf$node), size=1, shape=18)
  }
   p + geom_treescale(linesize=branchwidth, width = scalewidth, fontsize = 2)+
    geom_rootedge(.05, linewidth=branchwidth)+
    theme(legend.position = c(.15,.3),
          legend.background = element_rect(fill='transparent', color=NA),
          legend.box.background = element_rect(fill='transparent', color=NA),
          legend.text.align = 0)
}


plot_minimal_phylotree <- function(tree, branchwidth=.25, scalewidth=0.5,
                           anno.size=1.5, ladderize=T, support=70,
                           label=label, group="Family", align=T)
{
  ggtree(tree, ladderize = ladderize, size=branchwidth)+
    geom_tiplab(aes(subset=Family!="NODE"), size=anno.size, align = align, linesize = 0.25)+
    geom_tiplab(aes(subset=Family=="NODE"), color='red', size=anno.size, align=align, linesize = 0.25)+
    geom_nodelab(aes(label=as.numeric(label), subset = !is.na(as.numeric(label)) & as.numeric(label) < 90), 
                 size=anno.size, hjust=-.01) +
    geom_treescale(linesize=branchwidth, width=scalewidth, fontsize=2)+
    geom_rootedge(.05, linewidth=branchwidth)+
    theme(legend.position = c(.15,.3),
          panel.background=element_blank(),
          plot.background=element_blank(),
          legend.background = element_rect(fill='transparent', color=NA),
          legend.box.background = element_rect(fill='transparent', color=NA),
          legend.text.align = 0)
}

#' function to italicize part of legends
toexpr <- function(x, plain = NULL) {
  getfun <- function(x) {
    ifelse(is.element(x, plain), "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

#' function to reduce legend size of plots
addSmallLegend <- function(myPlot, pointSize = 1, textSize = 2, spaceLegend = 0.1) 
  {myPlot+
    guides(shape = "none",
           color = guide_legend(override.aes = list(size = pointSize-.5, shape=15))) +
    theme(legend.title = element_text(size = textSize, face="bold"), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"),
          legend.text.align = 0,
          panel.background=element_blank(),
          plot.background=element_blank(),
          legend.background = element_blank(),
          legend.box.background = element_blank())
}