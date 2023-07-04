#' @import umap
#' @import rgl
Umap.RGL.GeneSets <- function(Object, main = "", labels = "", plot.legend = F, spheres = F)
{
  clear3d()
  Object.umap <- umap(Object@Data.RR,n_components = 3, random_state = max(Object@Data[[1]]$cluster))
  umap.labels <- (Object@Data[[1]]$cluster)

  colors=RColorBrewer::brewer.pal(n = length(unique(Object@Data[[1]]$cluster)), name = "Set1")



  layout <- Object.umap$layout
  layout <- cbind(layout, Object@Data[[1]]$Pathways)
  layout <- as.data.frame(layout)

  colnames(layout) <- c("Umap_D1", "Umap_D2", "umap_D3", "labels")
  layout$cluster <- Object@Data[[1]]$cluster
  if(spheres == T)
  {
    plot.type <- "s"
  }else{
    plot.type <- "p"
  }
  plotids <- plot3d(layout[,1],layout[,2], layout[,3], main = main,
                    xlab = colnames(layout)[1],
                    ylab = colnames(layout)[2],
                    zlab = colnames(layout)[3],#size = 20,
                               type=plot.type, col=as.integer(labels.u))
  if(!labels[1] == "")
  {
    layout[!layout$labels %in% labels,"labels"] <- ""
    plotlabs <- text3d(layout[,1],layout[,2], layout[,3],  layout$labels,
                       type="p", col=as.integer(labels.u))

  }
  if(plot.legend[1] ==T)
  {
    plot.legend <- legend3d("topright", legend = unique((Object@Data[[1]]$cluster)), pch = 16, col = as.integer(labels.u), cex=1, inset=c(0.02))
  }

  rglwidget(elementId = "plot3drgl" )


}
