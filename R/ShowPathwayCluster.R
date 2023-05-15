#' ShowPathwayCluster
#'
#' Show the GO terms and descriptions of the generated cluster with PlotPathwayCluster.
#' @import AnnotationDbi
#' @import GO.db
#'
#' @param Object a pathway object
#' @param cluster the number of the cluster want to show information or "ALL" for all the clusters information.
#' @param nPathway minimum number of pathways per group, used in PlotPathwayCluster
#'
#' @return data.frame
#'
#' @export
#'
setGeneric(name="ShowPathwayCluster",
           def=function(Object, cluster = "ALL", nPathway = 4)
           {
             standardGeneric("ShowPathwayCluster")
           }
)


#' ShowPathwayCluster
#'
#' @param Object a pathway object
#' @param nPathway minimum number of pathways per group
#'
#' @return data.frame
#'
#' @examples

setMethod(f="ShowPathwayCluster",
          signature = "PathwayObject",
          definition = function(Object, cluster="ALL", nPathway = 4)
{
  message("Make sure you are using the same nPathway value as in PlotPathwayCluster.")

  tryCatch(
    {
      if (toupper(cluster) != "ALL") {
        as.numeric(cluster)
      }
    },
    error=function(e) {
      message("The argument included for cluster parameter is not valid. Instead showing all the clusters.")
      cluster <- "ALL"
    }
  )

  df <- Object@DataPathways.RR
  maxgo <- max(df)
  ordergo <- hclust(1 - as.dist(df))$order
  mat_cor <- cor(df[ordergo, ordergo])

  res <- obtainDefCluster(mat_cor)

  if (length(res[which(lapply(res, function(x) length(x))>=nPathway)]) == 0) {
    message(paste0("There is no groups with at least ", nPathway, " pathways. Using the minimum value, 2 pathways per group."))
    go_id_list <- res[which(lapply(res, function(x) length(x))>=2)]
  } else {
    go_id_list <- res[which(lapply(res, function(x) length(x))>=nPathway)]
  }

  if (((cluster > length(go_id_list) | cluster <= 0)) & toupper(cluster) != "ALL")  {
    paste0(message("There is no cluster ", cluster, ". Instead showing all the clusters."))
    cluster = "ALL"
  }

  GOres <- lapply(go_id_list, obtainGOdescription)

  if (toupper(cluster) == "ALL") {
    finalRes <- GOres
  } else {
    finalRes <- GOres[[cluster]]
  }

  return(finalRes)
}
)

