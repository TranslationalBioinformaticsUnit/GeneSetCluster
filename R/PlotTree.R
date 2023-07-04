#' PlotTree
#'
#' Plots a tree plot with the distances per cluster
#' @import ggplot2
#' @import ggnewscale
#' @import ggtree
#' @import GO.db
#'
#' @param Object a pathway object
#' @param fontsize a numeric with the fontsize for plotting nodenames
#' @param main is the plot title
#' @param nodenames boolean to add names to the nodes
#' @param wordcloud boolean to add wordcloud annotations of each cluster
#'
#' @return plot
#'
#' @export
#'
setGeneric(name="PlotTree",
           def=function(Object, fontsize = 3,
                        main="",
                        nodenames=TRUE,
                        wordcloud=TRUE)
           {
             standardGeneric("PlotTree")
           }
)

#' PlotTree
#'
#' @param Object a pathway object
#' @param fontsize a numeric with the fontsize for plotting rownames/colnames
#' @param main is the plot title
#' @param nodenames boolean to add names to the nodes
#' @param wordcloud boolean to add wordcloud annotations of each cluster
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
#'                                      method = "Hierarchical")
#'  PlotTree(Object = man.Great.Object3, fontsize = 3,
#'               main= "man.Great.Object3",
#'               nodenames = TRUE,
#'               wordcloud = TRUE)

setMethod(f="PlotTree",
          signature="PathwayObject",
          definition=function(Object, fontsize = 3,
                              main="",
                              nodenames = TRUE,
                              wordcloud = TRUE)
          {


            getterm<- function(goid){
              termid<-GOTERM[[goid]]
              if(is.null(termid)) {
                return("NA")
              }else{
                return(Term(termid))
              }
            }
            a<-Object@Data[[1]][,"Pathways"]

            results<-as.data.frame(unlist(lapply(a,getterm)))
            colnames(results)<-"Term"

            clus.x<-hclust(dist(t(Object@Data.RR)), method = "ward.D2")
            #metadata<-as.data.frame(Object@Data[[1]][,c("RR_name","Ratio")])
            #metadata<-as.data.frame(Object@Data[[1]][,c("RR_name","Ratio","Pathways")])
            metadata<-as.data.frame(c(Object@Data[[1]][,c("RR_name","Ratio")], results))

            if (nodenames==FALSE){
              ###without rownames
              p<-ggtree(clus.x) %<+% metadata + geom_tippoint(aes(size=Ratio)) + theme_tree2() + ggtitle(main) + theme(plot.title = element_text(hjust=0.5, face="bold"))
              offset<-0.05
              offset.update<-offset+70
            }else{
              #p<-ggtree(clus.x) %<+% metadata + geom_tippoint(aes(size=Ratio)) + geom_tiplab(aes(label=Pathways),size=fontsize, offset=10, align=TRUE, linesize=0.2) + theme_tree2() + ggtitle(main) + theme(plot.title = element_text(hjust=0.5, face="bold"))
              p<-ggtree(clus.x) %<+% metadata + geom_tippoint(aes(size=Ratio)) + geom_tiplab(aes(label=Term),size=fontsize, offset=10, align=TRUE, linesize=0.2) + theme_tree2() + ggtitle(main) + theme(plot.title = element_text(hjust=0.5, face="bold"))
              offset<-max(nchar(rownames(Object@Data.RR)))*11
              offset.update<-(offset/2)
              offset<-max(nchar(metadata$Term))*11

            }
            for (i in 1:dim(Object@plot$aka2)[2]){
              p<-gheatmap(p,Object@plot$aka2[i], offset=offset, colnames_position="top", colnames_offset_y=0.05, width=0.05, font.size=fontsize) + scale_fill_manual(values=Object@plot$aka3[[i]], name=names(Object@plot$aka3)[i])
              offset<-offset+offset.update
              p<-p + new_scale_fill()
            }

            if (wordcloud == TRUE) {
            for (i in unique(Object@Data[[1]][, "cluster"])) {
              list <- Object@Data[[1]][,"Pathways"][Object@plot$aka2$Cluster==i]
              plot <- wordcloud_generation(list)
              plot <- plot + annotation_custom(grob = textGrob(label = paste0("Cluster ", i), hjust = 0.5, vjust=-3,
                                                               gp = gpar(fontsize = 18, fontface = "bold")))
              if(i==1) {
                word_plot <- plot
              } else {
                word_plot <- big_plot/plot
              }
            }

            plot_out <- plot_grid(p, word_plot, nrow = 1, ncol = 2,   rel_widths = c(3,1))
            } else {
              plot_out <- p
            }

            return(plot_out)
          }
)
