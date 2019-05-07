#' HighlightGeneSets
#'
#' Adds a highlight score if the Gene-Set overlaps with a gene subset which is supplied by the user.
#'
#' @param Object A PathwayObject.
#' @param highligt.genes A vector with genes from the subset the user is interested in. e.g. a list of ROS genes.
#' @param name The name of the subset which will be added to the score calculated.
#'
#' @return a pathwayobject
#' @export
#'
#' @examples
#' IPA.files <- c(system.file("extdata",
#'                            "MM10.IPA.KO.uGvsMac.Canonical_pathways.xls",
#'                             package = "GeneSetCluster"),
#'              system.file("extdata",
#'                             "MM10.IPA.WT.uGvsMac.Canonical_pathways.xls",
#'                              package = "GeneSetCluster"),
#'              system.file("extdata",
#'                              "MM10.IPA.KO.uGvsMac.Functional_annotations.xls",
#'                              package = "GeneSetCluster"),
#'              system.file("extdata",
#'                              "MM10.IPA.WT.uGvsMac.Functional_annotations.xls",
#'                              package = "GeneSetCluster"))
#' canonical.files <- IPA.files[grep("Canonical", IPA.files)]
#'
#' IPA.object1 <- LoadGeneSets(file_location = canonical.files,
#'                          groupnames= c("KO", "WT"),
#'                          P.cutoff = 1.3,
#'                          Mol.cutoff = 5,
#'                          Source = "IPA",
#'                          type = "Canonical_Pathways",
#'                          structure = "SYMBOL",
#'                          seperator = ",")
#' IPA.object2 <- CombineGeneSets(Object = IPA.object1)
#'
#' IPA.object3 <- ClusterGeneSets(Object = IPA.object2,
#'                               clusters = 12,
#'                               method = "kmeans")
#' system.file("data", "Redox.genes.rda", package = "testdat")
#' IPA.object3.highlight <- HighlightGeneSets(Object = IPA.object3,
#'                                           highligt.genes = Redox.genes,
#'                                           name = "Ros")
#'

HighlightGeneSets <- function(Object, highligt.genes, name = "Ros")
{
  message("[=========================================================]")
  message("[<<<<<             HighlightGeneSets start            >>>>>>]")

  if(is.na(unique(Object@metadata[,"cluster.method"])))
  {
    message("Make sure youre object has been clustered by ClusterGeneSets")
    stop()
  }

  message("make sure that the structure of the highlight.genes is the same as the data")
  #Make sure that the highlight genes have the same structure as the object
  if(unique(as.character(Object@metadata[,"structure"])) == "SYMBOL")
  {
    message("transforming all highligt.genes to upper case, make sure this doesnt change the data")
    message(paste( "raw data has " , length(unique(highligt.genes))," highligt.genes", sep=""))
    highligt.genes <- toupper(highligt.genes)
    message(paste( "Transformed data has " , length(unique(highligt.genes))," highligt.genes", sep=""))
  }

  ################################################
  #-----------Seperate per cluster---------------#
  ################################################

  message("Calculating overlap between pathways and highlight genes")


  data <- vector()
  for(cl.i in unique(Object@Data[[1]]$cluster))
  {
    DF.cl <- Object@Data[[1]][Object@Data[[1]]$cluster == cl.i,]

    highligt.score <- vector()
    for(genes.cl.i in 1:nrow(DF.cl))
    {
      genes.x <- as.vector(strsplit2(as.character(DF.cl[genes.cl.i,]$Molecules), split=Object@metadata[1,"seperator"]))
      y <- sum(highligt.genes %in% genes.x)
      highligt.score[genes.cl.i] <- y/length(genes.x)
    }
    DF.cl$Highlight <- highligt.score
    DF.cl$Highlight.mean <- mean(highligt.score)
    data <- rbind(data, DF.cl)
  }
  Object@Data <- list(data)
  ####################################################
  #-----------sort clusters per highlight------------#
  ####################################################

  Object@Data <- list(Object@Data[[1]][order(Object@Data[[1]]$Highlight.mean, decreasing = T),])
  Object@Data.RR <- Object@Data.RR[order(Object@Data[[1]]$Highlight.mean, decreasing = T),order(Object@Data[[1]]$Highlight.mean, decreasing = T)]
  Object@plot$aka2 <- Object@plot$aka2[order(Object@Data[[1]]$Highlight.mean, decreasing = T),]
  Object@metadata[,"order.group"] <- "highlight"


  Object@metadata[,"highlight"] <- name

  #######################################################
  #-----------ADD info to plot Object@plot-------------#
  #######################################################

  #This will add a bar to the plot called highlight, where the cluster mean of the highlighted genes are displayed in an heatmap fashion, the higher the highlight the darker blue the color
  Cols.hightlight <- round(as.numeric(as.character(Object@Data[[1]]$Highlight.mean))*100, digits = 0)+1
  if(length(unique(Cols.hightlight)) == 1)
  {
    message("No highlights found in dataset, try different set of highlight genes")
    stop()
  }

  if(max(Cols.hightlight) > 75)
  {
    col.ramp.highlight <- colorRampPalette(c("darkblue","Blue","skyblue", "lightblue","grey","grey85","grey95", "white"))
  }
  if(max(Cols.hightlight) > 50 & max(Cols.hightlight) < 75)
  {
    col.ramp.highlight <- colorRampPalette(c("Blue","skyblue", "lightblue","grey","grey85","grey95", "white"))
  }
  if(max(Cols.hightlight) > 35 & max(Cols.hightlight) < 50)
  {
    col.ramp.highlight <- colorRampPalette(c("skyblue", "lightblue","grey","grey85","grey95", "white"))
  }else{
    col.ramp.highlight <- colorRampPalette(c("lightblue","grey","grey85","grey95", "white"))

  }
  Pal.highlight <- col.ramp.highlight(101)
  pal.hightlight <- Pal.highlight[Cols.hightlight]


  Object@plot$aka2$Highlight <- (Cols.hightlight-1)
  Pal.highlight <- unique(pal.hightlight)
  names(Pal.highlight) <- unique((Cols.hightlight-1))
  Object@plot$aka3$Highlight <- Pal.highlight
  ##################################
  #-----------Return---------------#
  ##################################


  message("-----------------------------------------------------------")
  message("[<<<<<             HighlightGeneSets END              >>>>>>]")
  message("[=========================================================]")
  return(Object)
}
