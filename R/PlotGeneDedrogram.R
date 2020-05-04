#' PlotGeneDendrogram
#' @param Object a pathway object
#' @param load.order Boolean show to use the previous order of clusters, though this could result in non-hierarchical ordering of the dendrogram
#'
#' @return plot
#' @export
#'
#' @examples
#'IPA.files <- c(system.file("extdata", "MM10.IPA.KO.uGvsMac.Canonical_pathways.xls", package = "GeneSetCluster"),
#'               system.file("extdata", "MM10.IPA.WT.uGvsMac.Canonical_pathways.xls", package = "GeneSetCluster"),
#'               system.file("extdata", "MM10.IPA.KO.uGvsMac.Functional_annotations.xls",package = "GeneSetCluster"),
#'               system.file("extdata", "MM10.IPA.WT.uGvsMac.Functional_annotations.xls",package = "GeneSetCluster"))
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
#'  IPA.object2 <- CombineGeneSets(Object = IPA.object1)
#'  IPA.object3 <- ClusterGeneSets(Object = IPA.object2,
#'                             clusters = 7,
#'                             method = "kmeans")
#'
#'  PlotGeneDendrogram(Object = IPA.object3, labels=F, RRmin = 0)
#'
#' @export
#'
setGeneric(name="PlotGeneDendrogram",
           def=function(Object,
                        load.order =F)
           {
             standardGeneric("PlotGeneDendrogram")
           }
)

#' PlotGeneDendrogram
#' @param Object a pathway object
#' @param load.order Boolean show to use the previous order of clusters, though this could result in non-hierarchical ordering of the dendrogram
#'
#' @return plot
#' @export
#'
#' @examples
#'
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
#'  IPA.object2 <- CombineGeneSets(Object = IPA.object1)
#'  IPA.object3 <- ClusterGeneSets(Object = IPA.object2,
#'                             clusters = 7,
#'                             method = "kmeans")
#'
#'  PlotGeneDendrogram(Object = IPA.object3, load.order=F)
setMethod(f="PlotGeneDendrogram",
          signature="PathwayObject",
          definition=function(Object,
                              load.order =F)
          {
            clus.x <- (hclust(dist(t(Object@Data.RR)), method = "ward.D2"))

            if(load.order == T)
            {
              warning("Dendrogram created  with a custom order, the tree might not make sense")
              clus.x$order <- Object@Data[[1]]$cluster
            }
            plot(clus.x)
          }

)
