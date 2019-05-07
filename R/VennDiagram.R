#' VennDiagram
#'
#' Plots a Venndiagram of overlap. Based on the limma functions.
#'
#' @import limma
#' @param n_groups The number of groups to compare.
#' @param Group1 A vector.
#' @param Group2 A vector.
#' @param Group3 A vector.
#' @param Group4 A vector.
#' @param names_groups A vector with a names for each group.
#' @param main Plot title.
#' @param legend Add a legend to the plot.
#' @param percentage Add percentages to each group. Default is False
#'
#' @return VennDiagram
#'
#' @examples
#' VennDiagram(n_groups = 2,
#'             Group1 =c("A","B","C","D","E"),
#'             Group2 =c("A","B","F","G","H"),
#'             names_groups = c("KO", "WT"),
#'             main = "Overlap",
#'             legend = FALSE,
#'             percentage = FALSE )
#'
#' @export
#'
VennDiagram <- function(n_groups=4,Group1,Group2,Group3,Group4,names_groups, main, legend = F, percentage =F)
{
  if(n_groups==2)
  {
    #Names_groups is the name you want in the plot
    idx.cutoff <- unique(c(Group1, Group2))
    idx.cutoff.mtx <- matrix(data = 0, nrow = length(idx.cutoff), ncol = 2)
    colnames(idx.cutoff.mtx) <- names_groups
    rownames(idx.cutoff.mtx) <- idx.cutoff
    idx.cutoff.mtx[idx.cutoff %in% Group1,names_groups[1]] <- 1
    idx.cutoff.mtx[idx.cutoff %in% Group2,names_groups[2]] <- 1
    if(percentage == T)
    {
      vc <- vennCounts(idx.cutoff.mtx)
      temp <- vc[,"Counts"]/ sum(vc[,"Counts"]) *100
      vc[,"Counts"] <- paste(vc[,"Counts"], " (",round(temp, digits = 2), "%)", sep="")

    }else{
      vc <- vennCounts(idx.cutoff.mtx)
    }
    vennDiagram(object = vc, main = main, cex.main=1.7,cex = 1)
    if(legend ==T){
      legend("topright",legend=c(paste(names_groups[1]," = ",length(Group1)," Terms",sep=""),
                                 paste(names_groups[2]," = ",length(Group2)," Terms",sep="")), cex=0.6)
    }
  }
  if(n_groups == 3)
  {
    #Names_groups is the name you want in the plot
    idx.cutoff <- unique(c(Group1, Group2, Group3))
    idx.cutoff.mtx <- matrix(data = 0, nrow = length(idx.cutoff), ncol = 3)
    colnames(idx.cutoff.mtx) <- names_groups
    rownames(idx.cutoff.mtx) <- idx.cutoff
    idx.cutoff.mtx[idx.cutoff %in% Group1,names_groups[1]] <- 1
    idx.cutoff.mtx[idx.cutoff %in% Group2,names_groups[2]] <- 1
    idx.cutoff.mtx[idx.cutoff %in% Group3,names_groups[3]] <- 1
    if(percentage == T)
    {
      vc <- vennCounts(idx.cutoff.mtx)
      temp <- vc[,"Counts"]/ sum(vc[,"Counts"]) *100
      vc[,"Counts"] <- paste(vc[,"Counts"], " (",round(temp, digits = 2), "%)", sep="")

    }else{
      vc <- vennCounts(idx.cutoff.mtx)
    }
    vennDiagram(object = vc, main = main, cex.main=1.7,cex = 1)
    if(legend ==T){
      legend("topright",legend=c(paste(names_groups[1]," = ",length(Group1)," Terms",sep=""),
                                 paste(names_groups[2]," = ",length(Group2)," Terms",sep=""),
                                 paste(names_groups[3]," = ",length(Group3)," Terms",sep="")), cex=0.6)
    }
  }
  if(n_groups == 4)
  {
    #Names_groups is the name you want in the plot
    idx.cutoff <- unique(c(Group1, Group2, Group3, Group4))
    idx.cutoff.mtx <- matrix(data = 0, nrow = length(idx.cutoff), ncol = 4)
    colnames(idx.cutoff.mtx) <- names_groups
    rownames(idx.cutoff.mtx) <- idx.cutoff
    idx.cutoff.mtx[idx.cutoff %in% Group1,names_groups[1]] <- 1
    idx.cutoff.mtx[idx.cutoff %in% Group2,names_groups[2]] <- 1
    idx.cutoff.mtx[idx.cutoff %in% Group3,names_groups[3]] <- 1
    idx.cutoff.mtx[idx.cutoff %in% Group4,names_groups[4]] <- 1
    if(percentage == T)
    {
      vc <- vennCounts(idx.cutoff.mtx)
      temp <- vc[,"Counts"]/ sum(vc[,"Counts"]) *100
      vc[,"Counts"] <- paste(vc[,"Counts"], " (",round(temp, digits = 2), "%)", sep="")

    }else{
      vc <- vennCounts(idx.cutoff.mtx)
    }
    vennDiagram(object = vc, main = main, cex.main=1.7,cex = 1)
    if(legend ==T){
      legend("topright",legend=c(paste(names_groups[1]," = ",length(Group1)," Terms",sep=""),
                                 paste(names_groups[2]," = ",length(Group2)," Terms",sep=""),
                                 paste(names_groups[3]," = ",length(Group3)," Terms",sep=""),
                                 paste(names_groups[4]," = ",length(Group4)," Terms",sep="")), cex=0.6)
    }
  }
}


