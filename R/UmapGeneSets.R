#' UmapGeneSets
#'
#' @import stats
#' @import methods
#' @import limma
#'
#' @importFrom stats dist hclust kmeans umap
#'
#' @param Object A PathwayObject
#' @param main Header for the plot
#' @param labels Labels of GeneSets to be highlighted
#' @export

UmapGeneSets <- function(Object, main = "", labels = "")
{
  Object.umap <- umap(Object@Data.RR)
  umap.labels <- (Object@Data[[1]]$cluster)
  
  colors=brewer.pal(n = length(unique(Object@Data[[1]]$cluster)), name = "Set1")
  names(colors) <- unique(Object@Data[[1]]$cluster)
  
  
  layout <- Object.umap$layout
  
  par(mar=c(5.2,4.1,4.1,2.1))
  
  xlim <- range(layout[,1])
  ylim <- range(layout[,2])
  
  plot(xlim, ylim, type="n")
  
  points(layout[,1], layout[,2], col=colors[as.integer(umap.labels)],
         cex=1, pch=16)
  if(!labels[1] =="")
  {
    labels <- labels[labels %in% colnames(Object@Data.RR)]
    idx <-  which(colnames(Object@Data.RR) %in% labels)
    
    text(layout[idx,1], layout[idx,2], labels)
  }
  
  labels.u <- unique(umap.labels)
  legend.pos <- "bottomright"
  legend.text <- as.character(labels.u)
  
  legend(legend.pos, legend=legend.text, inset=0.03,
         col=colors[as.integer(labels.u)],
         pch=16)
}
