#' PlotGeneSets
#'
#' Plots a heatmap with the distances per cluster
#' @import pheatmap
#' @import RColorBrewer
#' @importFrom grDevices colorRampPalette
#'
#' @param Object a pathway object
#' @param fontsize a numeric with the fontsize for plotting
#' @param legend  add a legend to the plot
#' @param annotation.mol should the genes from the genes set be added to the plot.
#' @param main is the plot title
#' @param RR.max is the maximum distance size to be added, it cutsoff the max to this. Usefull if the distance score gets very high.
#' @param cluster.order is a user defined order of the clusters in the heatmap.
#'
#' @return plot
#'
#' @export
#'
setGeneric(name="PlotGeneSets",
           def=function(Object, fontsize = 5,
                        legend = T,
                        annotation.mol=F,
                        main="",
                        RR.max = "",
                        cluster.order = "")
           {
             standardGeneric("PlotGeneSets")
           }
)

#' PlotGeneSets
#'
#' @param Object a PathwayObject
#' @param PathwayObject a PathwayObject
#' @param fontsize a numeric with the fontsize for plotting
#' @param legend  add a legend to the plot
#' @param annotation.mol should the genes from the genes set be added to the plot.
#' @param main is the plot title
#' @param RR.max is the maximum distance size to be added, it cutsoff the max to this. Usefull if the distance score gets very high.
#' @param cluster.order is a user defined order of the clusters in the heatmap.
#'
#' @return plot
#' @export
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
#'  PlotGeneSets(Object = man.Great.Object3, fontsize =5,
#'               legend = TRUE,
#'               annotation.mol=FALSE,
#'               main="man.Great.Object3",
#'               RR.max = 50)
setMethod(f="PlotGeneSets",
          signature="PathwayObject",
          definition=function(Object, fontsize = 5,
                              legend = T,
                              annotation.mol=F,
                              main="",
                              RR.max = "",
                              cluster.order = "")
          {
            Data.RR <- Object@Data.RR
            if(!RR.max == "")
            {
              RR.max <- as.numeric(as.character(RR.max))
              for(rows.i in 1:nrow(Data.RR))
              {
                idx <- Data.RR[rows.i,] > RR.max
                names(idx) <- NULL
                Data.RR[rows.i,idx] <- RR.max
              }
            }
            if(!cluster.order[1]== "")
            {
              if(!length(cluster.order) == max(as.numeric(as.character(Object@plot$aka2$Cluster))))
              {
                message("Make sure youre the order of supplied clusters is the same length as the number of clusters")

              }else{
                order.x <- vector()
                for(order.i in cluster.order)
                {
                  order.x <- c(order.x, which(order.i == as.numeric(as.character(Object@plot$aka2$Cluster))))
                }
                Object@plot$aka2 <- Object@plot$aka2[order.x,]

                Data.RR <- Data.RR[order.x, order.x]
              }


            }

            col.ramp <- colorRampPalette(c("white","grey","Tomato","red","red","red","red","red","red","red"))
            if(annotation.mol==T)
            {
              pheatmap(mat = Data.RR, color = col.ramp(n = 75),
                       main=main,
                       fontsize_row = fontsize,
                       fontsize_col = fontsize,
                       labels_row = Object@Data[[1]]$Molecules,
                       labels_col = rownames(Object@plot$aka2),
                       cluster_cols = F,#Object@plot$col.dd,
                       cluster_rows = F,#Object@plot$col.dd,
                       cellwidth = fontsize,
                       cellheight = fontsize,
                       legend = legend,
                       annotation_col = Object@plot$aka2,
                       annotation_row = Object@plot$aka2,
                       annotation_colors = Object@plot$aka3,
                       annotation_legend = legend)
            }

            if(annotation.mol=="Pathway")
            {
              pheatmap(mat = Data.RR, color = col.ramp(n = 75),
                       main=main,
                       fontsize_row = fontsize,
                       fontsize_col = fontsize,
                       labels_row = rownames(Object@plot$aka2),
                       labels_col = rep("",times=nrow(Data.RR)),
                       cluster_cols = F,#Object@plot$col.dd,
                       cluster_rows = F,#Object@plot$col.dd,
                       cellwidth = fontsize,
                       cellheight = fontsize,
                       legend = legend,
                       annotation_col = Object@plot$aka2,
                       annotation_row = Object@plot$aka2,
                       annotation_colors = Object@plot$aka3,
                       annotation_legend = legend)
            }
            if(annotation.mol ==F)
            {
              pheatmap(mat = Data.RR,
                       color = col.ramp(n = 75),
                       main=main,
                       labels_row = rep("",times=nrow(Data.RR)),
                       labels_col = rep("",times=ncol(Data.RR)),
                       cluster_cols = F,#Object@plot$col.dd,
                       cluster_rows = F,#Object@plot$col.dd,
                       annotation_col = Object@plot$aka2,
                       annotation_row = Object@plot$aka2, #border_color="blue",
                       annotation_colors = Object@plot$aka3,
                       annotation_legend = T)
            }
          }
)

