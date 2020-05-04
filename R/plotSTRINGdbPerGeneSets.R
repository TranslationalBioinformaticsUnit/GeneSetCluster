#' plotSTRINGdbPerGeneSets
#'
#' @import STRINGdb
#'
#' @importFrom stats dist hclust kmeans
#' @param StringObject The output from GetSTRINGdbPerGeneSets
#' @param plot.cluster The cluster for which the string image should be plotted.
#'
#' @export
#' @examples
#'
#'IPA.files <- c(system.file("extdata", "MM10.IPA.KO.uGvsMac.Canonical_pathways.xls", package = "GeneSetCluster"),
#'system.file("extdata", "MM10.IPA.WT.uGvsMac.Canonical_pathways.xls", package = "GeneSetCluster"),
#'system.file("extdata", "MM10.IPA.KO.uGvsMac.Functional_annotations.xls", package = "GeneSetCluster"),
#'system.file("extdata", "MM10.IPA.WT.uGvsMac.Functional_annotations.xls", package = "GeneSetCluster"))
#'canonical.files <- IPA.files[grep("Canonical", IPA.files)]
#'
#'IPA.object1 <- LoadGeneSets(file_location = canonical.files, #where are  the files
#'                            groupnames= c("KO", "WT"),#Names of the groups
#'                            P.cutoff = 1.3,
#'                            Mol.cutoff = 5,
#'                            Source = "IPA",
#'                            type = "Canonical_Pathways",
#'                            structure = "SYMBOL",
#'                            seperator = ",")
#'IPA.object2 <- CombineGeneSets(Object = IPA.object1)
#'IPA.object3 <- ClusterGeneSets(Object = IPA.object2,
#'                               clusters = 7,
#'                               method = "kmeans")
#'StringObject_unique <- GetSTRINGdbPerGeneSets(Object = IPA.object3,
#'                                              unique.per.cluster = TRUE ,
#'                                              plot.input = "All")
#' par(mfrow=c(3,3))
#' plotSTRINGdbPerGeneSets(StringObject = StringObject_all, plot.cluster = 1)
#' plotSTRINGdbPerGeneSets(StringObject = StringObject_all, plot.cluster = 2)
#' plotSTRINGdbPerGeneSets(StringObject = StringObject_all, plot.cluster = 3)
#' plotSTRINGdbPerGeneSets(StringObject = StringObject_all, plot.cluster = 4)
#' plotSTRINGdbPerGeneSets(StringObject = StringObject_all, plot.cluster = 5)
#' plotSTRINGdbPerGeneSets(StringObject = StringObject_all, plot.cluster = 6)
#' plotSTRINGdbPerGeneSets(StringObject = StringObject_all, plot.cluster = 7)



plotSTRINGdbPerGeneSets <- function(StringObject, plot.cluster=1)
{
  plot(1:(dim(StringObject[[plot.cluster]]$img)[2]), type = "n", xaxt = "n", yaxt = "n",
       xlab = "", ylab = "",
       ylim = c(1, dim(StringObject[[plot.cluster]]$img)[1]),
       xlim = c(1, (dim(StringObject[[plot.cluster]]$img)[2])), asp = 1)

  mtext(StringObject[[plot.cluster]]$input, cex = 0.7)

  mtext(StringObject[[plot.cluster]]$link, cex = 0.7,side = 1)

  rasterImage(StringObject[[plot.cluster]]$img, 1, 1,
              (dim(StringObject[[plot.cluster]]$img)[2]),
              dim(StringObject[[plot.cluster]]$img)[1])
}
