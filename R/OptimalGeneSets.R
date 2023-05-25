#' CombineGeneSets
#'
#' Calculate distances between the different experiments.
#' @import clustree
#' @import cluster
#' @import factoextra
#' @import GGally
#'
#' @param object A PathwayObject.
#' @param method Which method to determing optimal number of clusters. gap, elbow or silhouette.
#' @param max_cluster Max number of clusters to test
#' @param cluster_method kmeans or hcut. Which clustering method is used
#' @param main A string to be used as title in the plot
#' @param pathwayData Boolean to set Pathway or Cluster data
#'
#' @return a plot
#'
#' @examples
#'
#'require(GeneSetCluster)
#'IPA.files <- c(system.file("extdata", "MM10.IPA.KO.uGvsMac.Canonical_pathways.xls",
#'               package = "GeneSetCluster"),
#'               system.file("extdata", "MM10.IPA.WT.uGvsMac.Canonical_pathways.xls",
#'               package = "GeneSetCluster"),
#'               system.file("extdata", "MM10.IPA.KO.uGvsMac.Functional_annotations.xls",
#'               package = "GeneSetCluster"),
#'               system.file("extdata", "MM10.IPA.WT.uGvsMac.Functional_annotations.xls",
#'               package = "GeneSetCluster"))
#'canonical.files <- IPA.files[grep("Canonical", IPA.files)]
#'
#'IPA.object1 <- LoadGeneSets(file_location = canonical.files, #where are  the files
#'                            groupnames= c("KO", "WT"),
#'                            P.cutoff = 1.3,
#'                            Mol.cutoff = 5,
#'                            Source = "IPA",
#'                            type = "Canonical_Pathways",
#'                            structure = "SYMBOL",
#'                            seperator = ",")
#'IPA.object2 <- CombineGeneSets(Object = IPA.object1)
#'OptimalGeneSets(Object = IPA.object2, method = "elbow", max_cluster= 24,
#'                cluster_method = "kmeans", main= "Kmeans for 24 clusters")
#'
#'OptimalGeneSets(Object = IPA.object2, method = "elbow", max_cluster= 24,
#'                cluster_method = "kmeans", main= "Kmeans for 24 clusters", pathwayData = FALSE")
#'
#' @export
#'
#'
OptimalGeneSets <- function(Object, method, max_cluster, cluster_method, main, pathwayData = FALSE)
{

  message("[=========================================================]")
  message("[<<<<            OptimalGeneSets START               >>>>>]")
  message("-----------------------------------------------------------")
  if(method %in% c("Elbow", "elbow")){method <- "elbow"}
  if(method %in% c("Silhouette", "silhouette")){method <- "silhouette"}
  if(method %in% c("Gap", "gap")){method <- "gap"}
  if(cluster_method %in% c("HierArchical", "hierarchical", "Hierarchical", "hclust")){cluster_method <- "hcut"}

  message("Finding optimal number of clusters")
  message(paste("Clustering method = ", cluster_method))
  message(paste("Optimizing method = ", method))

  if (pathwayData == TRUE) {
    RRData <- Object@DataPathways.RR
  } else {
    RRData <- Object@Data.RR
  }

  ######################################
  ##---------elbow--------------------##
  ######################################

  if(method == "elbow")
  {
    plot.x <-   fviz_nbclust(RRData, get(cluster_method), method = "wss", k.max = max_cluster) +ggtitle(paste(main,": The Elbow method",sep=""))
  }
  ######################################
  ##---------gap----------------------##
  ######################################

  if(method == "gap")
  {
    gap_stat <- clusGap(RRData, get(cluster_method), K.max = max_cluster, B = 50)
    plot.x <- fviz_gap_stat(gap_stat) + ggtitle(paste(main,": Gap Statistics",sep=""))
  }

  ######################################
  ##---------silhoutte----------------##
  ######################################

  if(method == "silhouette")
  {
    plot.x <-  fviz_nbclust(x = RRData, get(cluster_method), method = "silhouette", k.max = max_cluster) + ggtitle(paste(main,": The Silhouette Plot",sep=""))
  }

  message("-----------------------------------------------------------")
  message("[<<<<<             OptimalGeneSets END              >>>>>>]")
  message("[=========================================================]")
  message("[You may want to process ClusterGeneSets next.            ]")

  return(plot(plot.x))

}
