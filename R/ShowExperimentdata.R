#' ShowExperimentdata
#'
#' Show Experiment data
#'
#' @param Object a PathwayObject
#'
#' @return a table with experiment data.
#' @export

setGeneric(name="ShowExperimentdata",
           def=function(Object)
           {
             standardGeneric("ShowExperimentdata")
           }
)

#' ShowExperimentdata
#'
#' @param Object a PathwayObject
#' @param PathwayObject a PathwayObject
#'
#' @return a table with experiment data.
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
#' ShowExperimentdata(Object =IPA.object1 )
setMethod(f="ShowExperimentdata",
          signature="PathwayObject",
          definition=function(Object)
          {
            return(Object@PData)
          }
)
