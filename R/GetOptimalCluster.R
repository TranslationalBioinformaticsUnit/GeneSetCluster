
#' GetOptimalCluster
#'
#' @param Object A PathwayObject
#' @param max_cluster A numeric with the number of maximal possible clusters
#' @param method The statistical method for determining the optimal number, currently only supports gap statistics
#' @param cluster_method kmeans of hclust
#'
#'
#' @return k number of clusters
#'
#' @examples
#' GetOptimalGeneSetsClusters(Object, max_cluster = 6)
#' @export
#'

GetOptimalCluster <- function(Object, cluster_method = "kmeans",max_cluster = 15, method ="gap")
{

if (method %in% c("Gap", "gap")) {
  method <- "gap"
}
if (cluster_method %in% c("HierArchical", "hierarchical", 
                          "Hierarchical", "hclust")) {
  cluster_method <- "hcut"
}

  gap_stat <- clusGap(Object@Data.RR, get(cluster_method), 
                      K.max = max_cluster, B = 50)
  
  gap <- gap_stat$Tab[, "gap"]
  se <- gap_stat$Tab[, "SE.sim"]
  decr <- diff(gap) <= 0
  maxSE = list(method = "firstSEmax", 
               SE.factor = 1)
  
  #Method
  ##########################################################################  
  .maxSE <- function (f, SE.f, method = c("firstSEmax", "Tibs2001SEmax", 
                                          "globalSEmax", "firstmax", "globalmax"), SE.factor = 1) 
  {
    method <- match.arg(method)
    stopifnot((K <- length(f)) >= 1, K == length(SE.f), SE.f >= 
                0, SE.factor >= 0)
    fSE <- SE.factor * SE.f
    switch(method, firstmax = {
      decr <- diff(f) <= 0
      if (any(decr)) which.max(decr) else K
    }, globalmax = {
      which.max(f)
    }, Tibs2001SEmax = {
      g.s <- f - fSE
      if (any(mp <- f[-K] >= g.s[-1])) which.max(mp) else K
    }, firstSEmax = {
      decr <- diff(f) <= 0
      nc <- if (any(decr)) which.max(decr) else K
      if (any(mp <- f[seq_len(nc - 1)] >= f[nc] - fSE[nc])) which(mp)[1] else nc
    }, globalSEmax = {
      nc <- which.max(f)
      if (any(mp <- f[seq_len(nc - 1)] >= f[nc] - fSE[nc])) which(mp)[1] else nc
    })
  }
  ##########################################################################  
  k <- .maxSE(gap, se, method = maxSE$method, SE.factor = maxSE$SE.factor)
  
  
  return(k)
  
}
