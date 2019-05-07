#' WriteGeneSets
#'
#' Exports the pathways into a csv file
#'
#' @param Object is a PathwayObject.
#' @param file_location where to write the files to
#' @param name name to be added to the files, what experiments are these
#' @param write what to write, either "Data", "RR" or "Both"
#'
#' @return written tables in the folder designated.
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
#' WriteGeneSets(Object= IPA.object3,
#'              file_location = "~/Project9/test/",
#'              name = "IPA_20181123", write = "Both")
#'
#'
setGeneric(name="WriteGeneSets",
           def=function(Object, file_location = "~/Project9/test/", name = "IPA_20181123", write = "Both")
           {
             standardGeneric("WriteGeneSets")
           }
)
#' WriteGeneSets
#'
#' Exports the pathways into a csv file
#'
#' @import utils
#' @param Object is a PathwayObject.
#' @param PathwayObject  a PathwayObject
#' @param file_location where to write the files to
#' @param name name to be added to the files, what experiments are these
#' @param write what to write, either "Data", "RR" or "Both"
#'
#' @return written tables in the folder designated.
#' @export
#'
#' @examples
#' Wx <- 1
#'
setMethod(f="WriteGeneSets",
          signature="PathwayObject",
          definition=function(Object, file_location = "~/Project9/test/", name = "IPA_20181123", write = "Both")
          {
            if(write=="Data")
            {
              write.table(Object@Data[[1]],
                          file= paste(file_location,name,"_Pathway",".csv",sep=""),
                          sep=";",
                          row.names = F)
            }
            if(write=="RR")
            {
              write.table(Object@Data.RR,
                          file= paste(file_location,name,"_RR",".csv",sep=""),
                          sep=";",
                          row.names = F)
            }
            if(write=="Both")
            {
              write.table(Object@Data[[1]],
                          file= paste(file_location,name,"_Pathway",".csv",sep=""),
                          sep=";",
                          row.names = F)
              write.table(Object@Data.RR,
                          file= paste(file_location,name,"_RR",".csv",sep=""),
                          sep=";",
                          row.names = F)
            }
          }
)
