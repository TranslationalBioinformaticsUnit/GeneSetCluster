#' MergeObjects
#'
#' Merge 2 different PathwayObjects
#'
#' @param Object.1 A PathwayObject.
#' @param Object.2 A PathwayObject.
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
#' Great.files.NoBCKGRND <- Great.files[!grepl("BCKGRND", Great.files)]
#'
#' Great.Object1 <- LoadGeneSets(file_location = Great.files.NoBCKGRND,
#'                                       groupnames= c("KO", "WT"),
#'                                       P.cutoff = 0.05,
#'                                       Mol.cutoff = 5,
#'                                       Source = "Great",
#'                                       Great.Background = FALSE,
#'                                       type = "Canonical_Pathways",
#'                                     topranks = 20,
#'                                    structure = "SYMBOL",
#'                                    Organism = "org.Mm.eg.db",
#'                                    seperator = ",")
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
#' GREAT.merged.object <- MergeObjects(Object.1 =Great.Object1,
#'                                    Object.2 = Great.bckgnrd.Object1)
#'
#' @export
#'
MergeObjects <- function(Object.1, Object.2)#combine multiple loadings into 1 Object
{
  message("[=========================================================]")
  message("[<<<<            MergeObjects START                 >>>>>]")
  message("-----------------------------------------------------------")

  if(!nrow(Object.1@Data.RR)==0 |!nrow(Object.2@Data.RR)==0)
  {
    message("Make sure youre objects have not yet been combined or clustered")
    stop()
  }
  ######################################
  ##---------Merge data-------------##
  ######################################

  message("merging data sets")


  ######################################
  ##---------Combine data-------------##
  ######################################
  Data <- do.call(c, list(Object.1@Data, Object.2@Data))


  PData <- rbind(Object.1@PData, Object.2@PData)
  rownames(PData) <-  paste("Experiment",1:nrow(PData), sep="_")

  metadata <- rbind(Object.1@metadata, Object.2@metadata)
  rownames(metadata) <-  paste("Experiment",1:nrow(metadata), sep="_")



  Object <- PathwayObject(Data = (Data),
                          PData = PData,
                          metadata = metadata,
                          Data.RR = data.frame(),
                          plot = list(aka2 = data.frame(),
                                      aka3 = data.frame()))
  ################################
  ##---------Warnings-----------##
  ################################

  if(length(unique(Object@metadata[,"type"])) > 1)
  {
    message("Warning, data types not the same")
  }
  if(length(unique(Object@metadata[,"organism"])) > 1)
  {
    message("Warning, Not the same organism")
  }

  message("-----------------------------------------------------------")
  message("[<<<<<               MergeObjects END              >>>>>>]")
  message("[=========================================================]")

  return(Object)
}
