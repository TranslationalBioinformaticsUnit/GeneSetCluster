#' PlotLD
#'
#' Plots a LD with the distances per cluster
#' @import ggplot2
#' @import ggnewscale
#' @import reshape2
#'
#' @param Object a pathway object
#' @param fontsize a numeric with the fontsize for plotting rownames/colnames
#' @param annotation.fontsize a numeric with the fontsize for plotting the annotation labels
#' @param main is the plot title
#' @param label.triangle label name for plotting the triangles
#'
#' @return plot
#'
#' @export
#'
setGeneric(name="PlotLD",
           def=function(Object, fontsize = 8,
                        annotation.fontsize = 3,
                        main="",
                        label.triangle = "")
           {
             standardGeneric("PlotLD")
           }
)

#' PlotLD
#'
#' @param Object a pathway object
#' @param fontsize a numeric with the fontsize for plotting rownames/colnames
#' @param annotation.fontsize a numeric with the fontsize for plotting the annotation labels
#' @param main is the plot title
#' @param label.triangle label name for plotting the triangles
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
#'  PlotLD(Object = man.Great.Object3, fontsize =8,
#'               annotation.fontsize = 3,
#'               main="man.Great.Object3",
#'               label.triangle = "Cluster")

setMethod(f="PlotLD",
          signature="PathwayObject",
          definition=function(Object, fontsize = 8,
                              annotation.fontsize = 3,
                              main="",
                              label.triangle = "")
  {
  
    if (!label.triangle=="") {
      #order base on the labels of the triangle
      order.x<-order(Object@plot$aka2[[label.triangle]])
      Data.RR<-Object@Data.RR[order.x,order.x]
      Object@plot$aka2<-Object@plot$aka2[order.x,]
    }else{
      Data.RR<-Object@Data.RR
    }
    
    
    #generate correlation matrix (upper triangle)
    cordata<-round(cor(Data.RR),2)
    cordata[lower.tri(cordata)] <- NA
    melted_cordata<-melt(cordata, na.rm=TRUE)
    
    baseplot<-ggplot(data=melted_cordata, aes(x=Var1, y=Var2, fill=value)) + geom_tile(color="white") + ggtitle(main) +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation") +
      theme(plot.title = element_text(hjust=0.5, face="bold"), panel.background=element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(size=fontsize, angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(size=fontsize))
    #theme(plot.title = element_text(hjust=0.5, face="bold"), panel.background=element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
        
    limit<-dim(Data.RR)[1]
    rectmin<-seq(from=0.5, to=(limit-0.5))
    rectmax<-seq(from=1.5, to=(limit+0.5))
    posmin=-1.5
    posmax=0
    textpos=-0.7
    
    #add annotation bar plots
    for (i in 1:dim(Object@plot$aka2)[2]){
      if(!colnames(Object@plot$aka2)[i]==label.triangle){
        df<-data.frame(value=factor(Object@plot$aka2[[i]]), start=rectmin, end=rectmax, posmin=rep(posmin,limit), posmax=rep(posmax,limit))
        baseplot<- baseplot + new_scale_fill()
        baseplot<- baseplot + geom_rect(data=df, inherit.aes = FALSE, aes(xmin=start, xmax=end, ymin=posmin, ymax=posmax, fill=value)) + geom_rect(data=df, inherit.aes = FALSE, aes(xmin=posmin, xmax=posmax, ymin=start, ymax=end, fill=value)) + scale_fill_manual(values = Object@plot$aka3[[i]]) + 
          guides(fill=guide_legend(title=colnames(Object@plot$aka2)[i])) +
          annotate("text", x=textpos, y=textpos, label= colnames(Object@plot$aka2)[i], size = annotation.fontsize, color="black",angle = 0, fontface = "bold")
        baseplot<- baseplot
        posmax<-posmax-1.8
        posmin<-posmin-1.8
        textpos<-textpos-1.5
      }
    }
    
    if (!label.triangle=="") {
      ###add triangle annotation on the label of interest
      for (i in unique(Object@plot$aka2[[label.triangle]])){
        minval<-min(which(Object@plot$aka2[label.triangle]==i))
        maxval<-max(which(Object@plot$aka2[label.triangle]==i))
        distance<-(maxval-minval-1)/2
        baseplot<- baseplot + geom_segment(aes_string(x=minval-0.5,xend=minval+0.5,y=minval-0.5,yend=minval-0.5), color="black", linewidth=0.7) +
          geom_segment(aes_string(x=minval+0.5,xend=maxval+0.5,y=minval-0.5,yend=maxval-0.5), color="black", linewidth=0.7) +
          geom_segment(aes_string(x=minval-0.5,xend=maxval+0.5,y=maxval+0.5,yend=maxval+0.5), color="black", linewidth=0.7) +
          geom_segment(aes_string(x=minval-0.5,xend=minval-0.5,y=minval-0.5,yend=maxval+0.5), color="black", linewidth=0.7) +
          geom_segment(aes_string(x=maxval+0.5,xend=maxval+0.5,y=maxval-0.5,yend=maxval+0.5), color="black", linewidth=0.7) +
          annotate("text", x=minval+distance+3, y=minval+distance, label= i, size = annotation.fontsize, color="black",angle = 0, fontface = "bold")
      }
    }
    baseplot
  
  }
)