#' BreakUpCluster
#'
#' @import stats
#' @import methods
#' @import limma
#'
#' @importFrom stats dist hclust kmeans
#'
#' @param Object A PathwayObject
#' @param breakup.cluster An integer with the number of the cluster to be broken up
#' @param sub.cluster The amount of subclusters to be generated
#' @return PathwayObject
#' @export
#'
#' @examples
#' Great.files <- c(system.file("extdata", "MM10.GREAT.KO.uGvsMac.bed.tsv",
#'                              package = "GeneSetCluster"),
#' system.file("extdata", "MM10.GREAT.KO.uGvsMac.bed_BCKGRND.tsv", package = "GeneSetCluster"),
#' system.file("extdata", "MM10.GREAT.WT.uGvsMac.bed.tsv", package = "GeneSetCluster"),
#' system.file("extdata", "MM10.GREAT.WT.uGvsMac.bed_BCKGRND.tsv", package = "GeneSetCluster"))
#' Great.files.bckgrnd <- Great.files[grepl("BCKGRND", Great.files)]
#'
#'
#' Great.bckgnrd.Object1 <- LoadGeneSets(file_location = Great.files.bckgrnd,
#'                                       groupnames= c("KO", "WT"),
#'                                       P.cutoff = 0.05,
#'                                       Mol.cutoff = 5,
#'                                       Source = "Great",
#'                                       Great.Background = TRUE,
#'                                       type = "Canonical_Pathways",
#'                                     topranks = 20,
#'                                    structure = "SYMBOL",
#'                                    Organism = "org.Mm.eg.db",
#'                                    seperator = ",")
#' man.Great.Object1 <- ManageGeneSets(Object = Great.bckgnrd.Object1,
#'                                    keep.type =c("Disease Ontology",
#'                                    "GO Biological Process" ),
#'                                    exclude.type="")
#' man.Great.Object2 <- CombineGeneSets(Object = man.Great.Object1)
#' man.Great.Object3 <- ClusterGeneSets(Object = man.Great.Object2,
#'                                      clusters = 5,
#'                                      method = "kmeans")
#'                                      
#' man.Great.Object3 <- (Object = man.Great.Object3, breakup.cluster = 3, sub.cluster=3)                                     
#'                                      

BreakUpCluster <- function(Object = Object, breakup.cluster = 6, sub.cluster=3)
{
  
  message("[=========================================================]")
  message("[<<<<            BreakUpCluster START                >>>>>]")
  message("-----------------------------------------------------------")
  
  require(RColorBrewer)
  
  
  warning("Currently only working for Kmeans and Hierarchical")
  if(is.na(Object@metadata$cluster.method[1]))
  {
    message("Make sure youre object has been combined by CombineGeneSets")
    stop()
  }  
  canonical.df <- Object@Data[[1]]
  
  
    order <- Object@metadata$order.group[1]
  method <- Object@metadata$cluster.method[1]
  
  
  #############################################
  ###----Kmeans Cluster--------------------###
  #############################################
  
  if(Object@metadata$cluster.method[1] == "kmeans")
  {
    Data.RR.clus <- Object@Data.RR[Object@Data[[1]]$cluster==breakup.cluster,Object@Data[[1]]$cluster==breakup.cluster]
    canonical.df$cluster[Object@Data[[1]]$cluster==breakup.cluster] <- paste0(breakup.cluster, ".",kmeans(x = Data.RR.clus, centers = sub.cluster)$cluster)
    
    Object@metadata[,"cluster.method"] <- paste0(Object@metadata[,"cluster.method"], "_breakup")
  }
  
  #############################################
  ###----Hierarchical Cluster--------------------###
  #############################################
  
  if(Object@metadata$cluster.method[1] == "hierarchical")
  {
    Data.RR.clus <- Object@Data.RR[Object@Data[[1]]$cluster==breakup.cluster,Object@Data[[1]]$cluster==breakup.cluster]
    temp.clx <- cutree(hclust(dist(t(Data.RR.clus)), method = "ward.D2"), k = sub.cluster)
    canonical.df$cluster[Object@Data[[1]]$cluster==breakup.cluster] <- paste0(breakup.cluster, ".",temp.clx)
    
    Object@metadata[,"cluster.method"] <- paste0(Object@metadata[,"cluster.method"], "_breakup")
  }
  Object@Data[[1]] <- canonical.df
  ############################################
  #---------Adding cluster info--------------#
  ############################################
  
  canonical.df$mean.RR.cl <- NA
  canonical.df$sum.RR.cl <- NA
  
  
  for(cl.i in unique(canonical.df$cluster))
  {
    DF.cl <- canonical.df[canonical.df$cluster == cl.i,]
    idx <- DF.cl$RR_name
    rr.x <- Object@Data.RR[idx, idx]
    rr.x <- data.frame(sapply(rr.x, function(x) as.numeric(as.character(x))))
    rr.x <- as.matrix(rr.x)
    DF.cl$mean.RR.cl <- round(mean(rr.x),digits = 4)
    DF.cl$sum.RR.cl <- round(sum(rr.x),digits = 4)
    
    canonical.df[canonical.df$cluster == cl.i,] <- DF.cl
    
    
  }
  
  ###############################
  #--------Order clusters-------#
  ###############################
  message("Ordering pathway clusters")
  if(order == "group")
  {
    Object@Data <- list(canonical.df[order(canonical.df$Group,canonical.df$cluster,  decreasing = F),])
    
    Object@metadata[,"order.group"] <- rep(order, times = nrow(Object@metadata))
    
  }
  if(order == "cluster")
  {
    Object@Data <- list(canonical.df[order(canonical.df$cluster,  decreasing = F),])
    
    Object@metadata[,"order.group"] <- rep(order, times = nrow(Object@metadata))
  }
  if(order == "mean.RR")
  {
    Object@Data <- list(canonical.df[order(canonical.df$mean.RR.cl,  decreasing = F),])
    
    Object@metadata[,"order.group"] <- rep(order, times = nrow(Object@metadata))
    
  }
  if(order == "Sum.RR")
  {
    Object@Data <- list(canonical.df[order(canonical.df$sum.RR.cl,  decreasing = F),])
    
    Object@metadata[,"order.group"] <- rep(order, times = nrow(Object@metadata))
    
  }
  Object@metadata[,"cluster.method"] <- rep(method, times = nrow(Object@metadata))
  Object@Data.RR <- Object@Data.RR[Object@Data[[1]]$RR_name,Object@Data[[1]]$RR_name]
  
  ################################################
  ####-----------Plotting Info----------------####
  ################################################
  
  
  pal.c <-  c(brewer.pal(n = 8, name ="Accent" ),
              brewer.pal(n = 8, name ="Dark2" ),
              brewer.pal(n = 8, name ="Set3"),
              brewer.pal(n = 8, name ="Set1"))
  
  if(Object@metadata$display[1] == "Expanded")
  {
    #Display is expanded meaning groups get marked seperatly
    if(length(unique(as.character(Object@Data[[1]]$Type))) > 1)
    {
      aka2 <- matrix(data = NA,
                     nrow = nrow((Object@Data.RR)),
                     ncol = 3 +nrow(Object@metadata))
      colnames(aka2) <- c("Group", "Cluster", "Type", unique(Object@metadata$Groups))
      aka2[,"Type"] <- as.character(Object@Data[[1]]$Type)
      pal.type <- pal.c[1:length(unique(aka2[,"Type"]))]
      names(pal.type) <- unique(aka2[,"Type"])
      
      for(groups.i in 1:nrow(Object@metadata))
      {
        aka2[,Object@metadata$Groups[groups.i]] <- Object@Data[[1]][,Object@metadata$Groups[groups.i]]
        
      }
      
    }else{
      aka2 <- matrix(data = NA, nrow = nrow((Object@Data.RR)), ncol = 2+ nrow(Object@metadata))
      colnames(aka2) <- c("Group", "Cluster", unique(Object@metadata$Groups))
      
      for(groups.i in 1:nrow(Object@metadata))
      {
        aka2[,Object@metadata$Groups[groups.i]] <- Object@Data[[1]][,Object@metadata$Groups[groups.i]]
        
      }
    }
  }else{
    if(length(unique(as.character(Object@Data[[1]]$Type))) > 1)
    {
      aka2 <- matrix(data = NA, nrow = nrow((Object@Data.RR)), ncol = 3)
      colnames(aka2) <- c("Group", "Cluster", "Type")
      aka2[,"Type"] <- as.character(Object@Data[[1]]$Type)
      pal.type <- pal.c[1:length(unique(aka2[,"Type"]))]
      names(pal.type) <- unique(aka2[,"Type"])
      
    }else{
      aka2 <- matrix(data = NA, nrow = nrow((Object@Data.RR)), ncol = 2)
      colnames(aka2) <- c("Group", "Cluster")
    }
  }
  rownames(aka2) <- Object@Data[[1]]$RR_name
  
  aka2[,"Group"] <- as.character(Object@Data[[1]]$Groups)
  aka2[,"Cluster"] <- as.character(Object@Data[[1]]$cluster)
  
  aka2 <- as.data.frame(aka2)
  if(length(unique(Object@Data[[1]]$cluster)) > 32)
  {
    message("Warning: number of cluster is larger than colours supported")
  }
  
  groups.col <- brewer.pal(n = 8, name ="Set2" )[1:length(unique(aka2[,"Group"]))]
  names(groups.col) <- unique(aka2[,"Group"])
  pal <- pal.c[1:length(unique(aka2[,"Cluster"]))]
  names(pal) <- unique(aka2[,"Cluster"])
  
  if(Object@metadata$display[1] == "Expanded")
  {
    
    if(length(unique(as.character(Object@Data[[1]]$Type))) > 1)
    {
      aka3 = list(Group = groups.col,
                  Cluster= pal,
                  Type = pal.type)
      
      for(groups.i in 1:nrow(Object@metadata))
      {
        vector.x <- c("black","red")
        names(vector.x) <- c(0,1)
        aka3[[Object@metadata$Groups[groups.i]]] <-   vector.x
      }
      
      names(aka3) <-  c("Group", "Cluster", "Type",  unique(Object@metadata$Groups))
      
    }else{
      aka3 = list(Group = groups.col,
                  Cluster= pal)
      
      for(groups.i in 1:nrow(Object@metadata))
      {
        vector.x <- c("black","red")
        names(vector.x) <- c(0,1)
        aka3[[Object@metadata$Groups[groups.i]]] <-   vector.x
      }
      
      names(aka3) <-  c("Group", "Cluster",  unique(Object@metadata$Groups))
    }
  }else{
    
    if(length(unique(as.character(Object@Data[[1]]$Type))) > 1)
    {
      aka3 = list(Group = groups.col,
                  Cluster= pal,
                  Type = pal.type)
      names(aka3) <-  c("Group", "Cluster", "Type")
      
    }else{
      aka3 = list(Group = groups.col,
                  Cluster= pal)
      names(aka3) <-  c("Group", "Cluster")
    }
  }
  
  ##################################
  #-----------Return---------------#
  ##################################
  
  Object@plot$aka2 = aka2
  Object@plot$aka3 = aka3
  
  
  message("-----------------------------------------------------------")
  message("[<<<<<             BreakUpCluster END               >>>>>>]")
  message("[=========================================================]")
  message("[You may want to process HighlightGeneSets next.           ]")
  message("[You may want to plot the results using PlotGeneSets next. ]")
  
  
  return(Object)
}



