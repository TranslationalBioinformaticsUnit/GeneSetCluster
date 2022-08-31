#' BreakUpCluster
#'
#' @import stats
#' @import methods
#' @import limma
#'
#' @importFrom stats dist hclust kmeans
#'
#' @param Object A PathwayObject
#' @return PathwayObject
#' @export
MeanClusterDifferences <- function(Object)
{
  message("[=========================================================]")
  message("[<<<<      MeanClusterDifferences START              >>>>>]")
  message("-----------------------------------------------------------")
  
  RR.Data <- Object@Data.RR
  Ncluster <- max(Object@Data[[1]]$cluster)
  
  Cluster.Combo <- t(combn(1:Ncluster, 2))
  Cluster.Combo <- cbind(Cluster.Combo,
                         matrix(NA, nrow = nrow(Cluster.Combo), ncol = 3))
 colnames(Cluster.Combo) <- c("Cluster_A", "Cluster_B",
                              "Cluster_A_Mean", "Cluster_B_Mean",
                              "Diff_A-B")

  for(cl.i in 1:nrow(Cluster.Combo))
  {
    idx_A <- Object@Data[[1]][Object@Data[[1]]$cluster == Cluster.Combo[cl.i,1], "RR_name"]
    idx_B <- Object@Data[[1]][Object@Data[[1]]$cluster == Cluster.Combo[cl.i,2], "RR_name"]
    
    Cluster.Combo[cl.i,"Cluster_A_Mean"] <- mean(as.numeric(unlist(RR.Data[idx_A,idx_A])))
    Cluster.Combo[cl.i,"Cluster_B_Mean"] <- mean(as.numeric(unlist(RR.Data[idx_B,idx_B])))
    
    diffs.vc <- vector()
    for(A.i in 1:length(idx_A))
    {
      for(B.i in 1:length(idx_B))
      {
        diffs.vc <- c(diffs.vc,
                      RR.Data[c(idx_A[A.i]),c(idx_B[B.i])])
      }
    }
    
  
    Cluster.Combo[cl.i,"Diff_A-B"] <- mean(as.numeric(diffs.vc))
    
  }
 
 Cluster.Combo <- as.data.frame(Cluster.Combo)
 Cluster.Combo$Ratio_A_Diff <-  Cluster.Combo$`Diff_A-B` / Cluster.Combo$Cluster_A_Mean
 Cluster.Combo$Ratio_B_Diff <- Cluster.Combo$`Diff_A-B` / Cluster.Combo$Cluster_B_Mean
 
 message("-----------------------------------------------------------")
 message("[<<<<<      MeanClusterDifferences END              >>>>>>]")
 message("[=========================================================]")
 return(Cluster.Combo)
}
