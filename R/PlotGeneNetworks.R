#' PlotGeneNetworks
#'
#' Plots a heatmap with the distances per cluster
#' @import factoextra
#' @import GGally
#' @import network
#'
#' @param Object a pathway object
#' @param labels Boolean show labels or not, or a vector with labels to be shown
#' @param RRmin filter for different RR cutoff for what is an edge
#'
#' @return plot
#'
#' @examples
#' @export
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
#'  PlotGeneNetworks(Object = IPA.object3, labels=F, RRmin = 0)
#'
#' @export
#'
setGeneric(name="PlotGeneNetworks",
           def=function(Object,
                        labels =F,
                        RRmin = 10)
           {
             standardGeneric("PlotGeneNetworks")
           }
)

#' PlotGeneNetworks
#' @import factoextra
#' @import GGally
#'
#' @param Object a pathway object
#' @param labels Boolean show labels or not, or a vector with labels to be shown
#' @param RRmin filter for different RR cutoff for what is an edge.
#'
#' @return plot
#' @export
#'
#' @examples
#' @export
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
#'  PlotGeneNetworks(Object = IPA.object3, labels=F, RRmin = 0)
setMethod(f="PlotGeneNetworks",
          signature="PathwayObject",
          definition=function(Object,
                              labels =F,
                              RRmin = 10)
          {
            RR.data <- as.matrix(Object@Data.RR)
            rownames(RR.data) <- paste("Cluster", Object@Data[[1]]$cluster,sep="")
            colnames(RR.data) <- paste("Cluster", Object@Data[[1]]$cluster,sep="")

            for(i in 1:nrow(RR.data))
            {
              RR.data[i,RR.data[i,] <= RRmin] <- 0
            }

            RR.net = network(RR.data,  directed = TRUE)

            # Cluster affiliation
            x <- as.factor(paste("Cluster", Object@Data[[1]]$cluster,sep=""))
            RR.net %v% "Cluster" = as.character(x)

            # color palette
            y <- Object@plot$aka3$Cluster
            names(y) <- paste("Cluster", names(y),sep="")


            if(sum(labels %in% T) == 1)
            {
              labels = Object@Data[[1]]$Pathways
            }
            if(length(labels) > 1 & length(labels)<= ncol(RR.data))
            {
              labels.2 <- rep("", times =ncol(RR.data))
              labels.2[Object@Data[[1]]$Pathways %in% labels] <- labels
              labels <- labels.2
            }


            # network plot
            ggnet2(RR.net, color = "Cluster", palette = y, alpha = 0.75, size = 4, edge.alpha = 0.5, label = labels)
          }

)
