#' PlotClusterTissueExpression
#'
#' Plots the tissue expression of one or more clusters. Tissue expression is ordered based on the first cluster entered.
#'
#' @import dplyr
#' @import patchwork
#' @import ggplot2
#'
#' @param Object a PathwayObject
#' @param lclusters numeric vector giving which cluster(s) plot
#' @param all boolean to add all tissues
#'
#' @return dataframe of pathwayobjects with tissue expression per cluster
#' @export
#'

setGeneric(name="PlotClusterTissueExpression",
           def=function(Object, lclusters, all=FALSE)
           {
             standardGeneric("PlotClusterTissueExpression")
           }
)

#' PlotClusterTissueExpression
#'
#' @param Object a pathway object
#' @param lclusters numeric vector giving which cluster(s) plot
#' @param all boolean to add all tissues
#'
#' @return plot
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
#' man.Great.Object4 <- TissueExpressionPerGeneSet(man.Great.Object3)
#' PlotClusterTissueExpression(man.Great.Object4,
#'                             lclusters = 2,
#'                             all = F)
#' PlotClusterTissueExpression(man.Great.Object4,
#'                             lclusters = c(2, 4, 5),
#'                             all = F)
#' PlotClusterTissueExpression(man.Great.Object4,
#'                             lclusters = c(2, 4, 5, 1, 3),
#'                             all = F)

setMethod(f="PlotClusterTissueExpression",
          signature = "PathwayObject",
          definition = function(Object,
                                lclusters,
                                all = F)
{

  if (is.null(Object@dfTissue) | anyNA(Object@dfTissue)) {
    stop("First you have to run TissueExpressionPerGeneSet to obtain tissue expression information.")
  }

  if (max(lclusters) > ncol(Object@dfTissue))
  {
    message(paste("Please select generated cluster(s) from 1 to "), ncol(Object@dfTissue))
    stop()
  }

  if (length(lclusters) == 0)
  {
    stop("Please include at least one cluster in the lclusters argument.")
  }

  #set maximum value
  max.value <- max(Object@dfTissue[, lclusters])
  lplots <- list()

  #select showing top tissues depending on showing number of clusters
  if (length(lclusters) < 4){
    top <- 20
  } else if (length(lclusters) < 6) {
    top <- 15
  } else {
    top <- 10
  }

  message("Note that the most expressed tissues are reordered according to the values of the first cluster.")
  message(paste("Therefore the most expressed top ", top, " tissues of cluster", lclusters[[1]], "are shown."), sep="")

  #generate plot for every cluster
  for (cluster in lclusters){
    data = as.data.frame(Object@dfTissue[, cluster])
    data$tissue = rownames(data)
    colnames(data)=c("value", "tissue")
    plot_name <- paste("Cluster", cluster)

    # reorder values according to the most exressed values in the first cluster
    if (lclusters[[1]]==cluster) {
      data = data[order(data$value, decreasing = T),]

      if (all == FALSE){
        data = data[1:top,] # top
        keep <- data$tissue

      }else{
        keep <- data$tissue # all
      }
    }

    data_plot <- data[rev(keep),] %>%
      mutate(tissue=factor(tissue, tissue))

    if (lclusters[[length(lclusters)]]==cluster)
    {
      p <- ggplot(data_plot,  aes(x=tissue, y=value) ) +
        geom_segment(aes(x=tissue ,xend=tissue, y=0, yend=value), color="black") +
        geom_point(size=3, color="#69b3a2") +
        ggtitle(plot_name) + ylab("median gene expression") + ylim(0, max.value) +
        coord_flip()

    }else{
      p <- ggplot(data_plot,  aes(x=tissue, y=value) ) +
        geom_segment(aes(x=tissue ,xend=tissue, y=0, yend=value), color="black") +
        geom_point(size=3, color="#69b3a2") +
        ggtitle(plot_name) +  ylim(0, max.value) +
        theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())  +
        coord_flip()
    }

    if (lclusters[[1]]==cluster)
    {
      p2 <- p

    }else{
      p2 <- p2/p
    }
  }

  return(p2)
}
)
