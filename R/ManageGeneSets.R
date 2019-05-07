#' ManageGeneSets
#'
#' Filter unwanted groups from the PathwayObject
#'
#' @param Object A PathwayObject.
#' @param exclude.type A vector with strings of groups meant to be removed from the PathwayObject while keeping all others.
#' @param keep.type A vector with strings of groups meant to be kept from the PathwayObject removing all others.
#'
#' @return a pathwayobject
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
#'
#' @export
#'
ManageGeneSets <- function(Object, exclude.type="", keep.type="")
{

  message("[=========================================================]")
  message("[<<<<<             ManageGeneSets start             >>>>>>]")

  if(nrow(Object@Data.RR) > 0)
  {
    message("Make sure youre objects have not yet been combined or clustered")
    stop()
  }

  #this function is meant to either exclude types or include types.
  #Especially for great, which has 18 potential types, its sometimes better to just look at 1 output
  if(!sum(exclude.type %in% "") >= 1)
  {
    message("excluding the following pathways types:" )
    for(types.i in 1:length(exclude.type))
    {
      message((exclude.type[types.i]))
    }
    if(sum(Object@metadata$type %in% exclude.type) < 1)
    {
      message("exclude.type not found")
      stop("try with other exclude.type")
    }
    exclude.meta <- Object@metadata[Object@metadata$type %in% exclude.type,]
    exclude.PData <- Object@PData[rownames(Object@PData) %in% rownames(Object@metadata[Object@metadata$type %in% exclude.type,]),]
    idx <- rep(F, length(names(Object@Data)))
    for(i in 1:length(exclude.type))
    {
      idx2 <- grepl(exclude.type[i], names(Object@Data))
      idx <- idx | idx2
    }
    exclude.Data <- Object@Data[idx]

    Object <- list(Data = exclude.Data,
                   PData = exclude.PData,
                   metadata = exclude.meta)
  }
  if(!sum(keep.type %in% "") >= 1)#so if keep.type is occupied
  {
    message("Keeping the following pathways types:")
    for(types.i in 1:length(keep.type))
    {
      message((keep.type[types.i]))
    }
    if(sum(Object@metadata$type %in% keep.type) < 1)
    {
      message("keep.type not found")
      stop()
    }
    exclude.meta <- Object@metadata[Object@metadata$type %in% keep.type,]
    exclude.PData <- Object@PData[rownames(Object@PData) %in% rownames(Object@metadata[Object@metadata$type %in% keep.type,]),]
    idx <- rep(F, length(names(Object@Data)))
    for(i in 1:length(keep.type))
    {
      idx2 <- grepl(keep.type[i], names(Object@Data))
      idx <- idx | idx2
    }
    exclude.Data <- Object@Data[idx]

    Object@Data = (exclude.Data)
    Object@PData = exclude.PData
    Object@metadata = exclude.meta

  }
  message("-----------------------------------------------------------")
  message("[<<<<<             ManageGeneSets END               >>>>>>]")
  message("[=========================================================]")

  return(Object)
}
