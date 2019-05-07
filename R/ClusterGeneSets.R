
#' ClusterGeneSets
#'
#' @import mclust
#' @import stats
#' @import methods
#' @importFrom stats dist hclust kmeans
#' @param Object A PathwayObject
#' @param clusters A numeric with the number of clusters required
#' @param method The cluster method specified
#' @param order How should the data be ordered, by group or by cluster
#'
#'
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
#'
#'

ClusterGeneSets <- function(Object, clusters = 5, method = "kmeans", order = "group")
{
  message("[=========================================================]")
  message("[<<<<            ClusterGeneSets START               >>>>>]")
  message("-----------------------------------------------------------")

  if(!nrow(Object@Data.RR)>=1)
  {
    message("Make sure youre object has been combines by CombineGeneSets")
    stop()
  }

  #Currently supported ways of clustering
  #Kmeans  #kmeans_group  #mclust #mclust_group #Hierarchical #Hierarchical_group
  #Correct spelling mistakes#
  if(method %in% c("kmeans", "Kmeans")){method <- "kmeans"}
  if(method %in% c("kmeans_group","Kmeans_group","Kmeans_Group", "kmeans group","Kmeans group")){method <- "kmeans_group"}
  #
  if(method %in% c("MClust", "Mclust", "mclust")){method <- "mclust"}
  if(method %in% c("Mclust_group","MClust_group", "mclust group","MClust group", "mclust_group")){method <- "mclust_group"}
  #
  if(method %in% c("HierArchical", "hierarchical", "Hierarchical")){method <- "Hierarchical"}
  if(method %in% c("Hierarchical_group","hierarchical_group", "Hierarchical group","hierarchical group", "hierarchical_group")){method <- "Hierarchical_group"}

  #####################################################
  if(clusters > nrow(Object@Data[[1]]))
  {
    message("Too many cluster inputted, reducing it to the number of datapoints")
    clusters <- nrow(Object@Data[[1]])
  }

  if(sum(method %in% c("kmeans", "kmeans_group", "mclust","mclust_group", "Hierarchical", "Hierarchical_group")) < 1)
  {
    message("Method is not a proper method of clustering")
    stop()
  }
  if(length(method) > 1)
  {
    message("Too many methods supplied, pick between Kmeans, kmeans_group, mclust, mclust_group")
    stop()
  }

  ####################################
  #------------Combine DF------------#
  ####################################

  canonical.df <-Object@Data[[1]]

  ####################################################
  #------------Kmeans of clusters--------------------#
  ####################################################
  if(method == "kmeans")
  {
    message("Running kmeans")

    if(length(clusters) > 1){message("Too many clusters supplied, Taking first cluster")}
    canonical.df$cluster <- kmeans(x = Object@Data.RR, centers = clusters)$cluster

  }


  ##############################################################
  #------------Kmeans of clusters per group--------------------#
  ##############################################################
  if(method == "kmeans_group")
  {
    print("Running kmeans per group seperatly")
    if(!length(clusters) == nrow(Object@Pdata))
    {
      message("Wrong number of clusters supplied")
      stop()
    }

    temp.cl2 <- vector()
    for(i in 1:nrow(Object@Pdata))
    {

      idx.ls <- Object@Data[[1]][as.character(Object@Data[[1]]$Groups) == as.character(Object@Pdata[i,"Groupnames"]),"RR_name"]

      rr.temp <- Object@Data.RR[idx.ls,idx.ls]

      temp.cl1 <-  kmeans(x = rr.temp, centers = clusters[i])$cluster
      if(i >1){temp.cl1 <- temp.cl1 + max(temp.cl2)}
      temp.cl2 <- c(temp.cl2, temp.cl1)
    }
    canonical.df$cluster <- temp.cl2


  }

  #############################################
  #----------mclust of clusters---------------#
  #############################################


  if(method == "mclust")
  {
    message("Running mclust")
    if(length(clusters) > 1)
    {
      message("Too many clusters supplied, Taking first cluster")
    }else{
      message(paste("maximum number of inputted clusters =", clusters))
    }

    mclust.Bic <- mclustBIC(Object@Data.RR,seq(from=2,to=clusters,by=1),
                            modelNames=c("VII", "EEI", "EVI", "VEI", "VVI", "EII"))

    mclust.BicNsummary <- summary(mclust.Bic, Object@Data.RR)
    canonical.df$cluster <- mclust.BicNsummary$classification
    print(paste("Optimal number of clusters =", max(mclust.BicNsummary$classification)))
  }



  ##############################################################
  #------------MClust of clusters per group--------------------#
  ##############################################################
  if(method == "mclust_group")
  {
    message("Running mclust on every group seperatly")
    if(!length(clusters) == nrow(Object@Pdata))
    {
      message("Wrong number of clusters supplied")
      stop()
    }
    optimal.clus <- vector()
    temp.cl2 <- vector()
    for(i in 1:nrow(Object@Pdata))
    {

      idx.ls <- Object@Data[[1]][as.character(Object@Data[[1]]$Group) == as.character(Object@Pdata[i,"Groupnames"]),"RR_name"]

      rr.temp <- Object@Data.RR[idx.ls,idx.ls]
      temp.BIC <- mclustBIC(rr.temp,seq(from=2,to=clusters[i],by=1),
                            modelNames=c("VII", "EEI", "EVI", "VEI", "VVI", "EII"))

      temp.BIC.sum <- summary(temp.BIC, rr.temp)
      temp.cl1 <- temp.BIC.sum$classification
      optimal.clus[i] <- max(temp.cl1)
      if(i >1){temp.cl1 <- temp.cl1 + max(temp.cl2)}
      temp.cl2 <- c(temp.cl2, temp.cl1)
    }
    canonical.df$cluster <- temp.cl2
    print(paste("Optimal number of clusters =",optimal.clus))

  }
  #########################################################
  #------------Hierarchical clustering--------------------#
  #########################################################
  if(method == "Hierarchical")
  {
    message("Running Hierarchical")

    if(length(clusters) > 1){message("Too many clusters supplied, Taking first cluster")}
    canonical.df$cluster <- hclust(dist(t(Object@Data.RR)), method = "ward.D2")$order

  }

  #########################################################
  #------------Hierarchical clustering--------------------#
  #########################################################
  if(method == "Hierarchical_group")
  {
    message("Running Hierarchical per group")

    if(length(clusters) > 1){message("Too many clusters supplied, Taking first cluster")}
    temp.cl2 <- vector()
    for(i in 1:nrow(Object@Pdata))
    {

      idx.ls <- Object@Data[[1]][as.character(Object@Data[[1]]$Group) == as.character(Object@Pdata[i,"Groupnames"]),"RR_name"]

      rr.temp <- Object@Data.RR[idx.ls,idx.ls]

      temp.cl1 <- hclust(dist(t(rr.temp)), method = "ward.D2")$order
      if(i >1){temp.cl1 <- temp.cl1 + max(temp.cl2)}
      temp.cl2 <- c(temp.cl2, temp.cl1)
    }
    canonical.df$cluster <- temp.cl2

  }
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

  ##########################
  #--------plot Info-------#
  ##########################

  pal.c <-  c(brewer.pal(n = 8, name ="Accent" ),
              brewer.pal(n = 8, name ="Dark2" ),
              brewer.pal(n = 8, name ="Set3"),
              brewer.pal(n = 8, name ="Set1"))


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

  ##################################
  #-----------Return---------------#
  ##################################

  Object@plot$aka2 = aka2
  Object@plot$aka3 = aka3



  message("-----------------------------------------------------------")
  message("[<<<<<             ClusterGeneSets END              >>>>>>]")
  message("[=========================================================]")
  message("[You may want to process HighlightGeneSets next.            ]")
  message("[You may want to plot the results using PlotGeneSets next. ]")


  return(Object)
}
