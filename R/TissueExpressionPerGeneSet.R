#' TissueExpressionPerGeneSet
#'
#' Extracts every gen of every cluster from the object. Using the median gene expression from the GTEx database build a data.frame with every tissue per cluster.
#'
#' @import jsonlite
#' @import doParallel
#' @import parallel
#' @import httr
#' @import reshape2
#' @import dplyr
#' @import utils
#'
#' @param Object a PathwayObject
#' @param PathwayObject a PathwayObject
#' @param threads number of cores. Default the maximum of the computer.
#'
#' @return dataframe of pathwayobjects with tissue expression per cluster
#'

setGeneric(name="TissueExpressionPerGeneSet",
           def=function(Object, threads=1)
           {
             standardGeneric("TissueExpressionPerGeneSet")
           }
)

#' TissueExpressionPerGeneSet
#'
#' @param Object a PathwayObject
#' @param PathwayObject a PathwayObject
#' @param threads number of cores. Default the maximum of the computer.
#'
#' @return dataframe of pathwayobjects with tissue expression per cluster
#' @export
#'
#' @examples
#' #' IPA.files <- c(system.file("extdata",
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
#' IPA.object3 <- ClusterGeneSets(Object = IPA.object2,
#'                               clusters = 7,
#'                               method = "kmeans")
#' IPA.object4 <- TissueExpressionPerGeneSet(Object = IPA.object3, threads = 8)
setMethod(f="TissueExpressionPerGeneSet",
          signature="PathwayObject",
          definition=function(Object, threads=1)
          {
          message("[=========================================================]")
          message("[<<<<         ObtainTissueExpression START           >>>>>]")
          message("-----------------------------------------------------------")


          if (is.na(threads)) {
            message(paste("Number of cores not specified, using "), detectCores(), " cores.")
            threads <- detectCores()

          } else if (threads > detectCores()) {
            message(paste("Introduced a higher number of cores than detected. Using the maximum of cores: "), detectCores())
            threads <- detectCores()
          }

          message("Please make sure that you are using GENE SYMBOLs.")

          data("dic", package = "GeneSetCluster")

          mol.unique.df <- GenesPerGeneSet(Object)
          genes <- rownames(mol.unique.df)

          if (substr(genes[1], 1, 4) == "ENSG") {
            dic.selected <- dic[which(dic$GTEx.median.TPM.Name %in% toupper(genes)),]
            genes.selected <- dic.selected$GTEx.median.TPM.Name
          } else {
          dic.selected <- dic[which(dic$GTEx.median.TPM.Description %in% toupper(genes)),]
          genes.selected <- dic.selected$GTEx.median.TPM.Name
          }

          if (length(genes.selected) == 0) {
            stop("The genes used to create the object are not correct. Please make sure they are introduced either as GENE SYMBOL or ENSEMBL ID.")
          }

          #GTEx API REST query of median Gene Expression: database -> gtex_v8
          message(paste("Performing GTEx API REST queries for ", length(genes.selected), " genes. Using ", threads, " cores.", sep=""))
          urls <- paste("https://gtexportal.org/rest/v1/expression/medianGeneExpression?datasetId=gtex_v8&gencodeId=", genes.selected, sep="")

          #parallelization settings
          cl <- makeCluster(threads)
          clusterEvalQ(cl, library(httr))
          results <- parLapply(cl, urls, httr::GET)
          stopCluster(cl)

          #extract data from json build the data frame
          results_df <- lapply(results, function(x) as.data.frame(cbind(fromJSON(rawToChar(x$content))$medianGeneExpression$gencodeId,
                                                                        fromJSON(rawToChar(x$content))$medianGeneExpression$geneSymbol,
                                                                        fromJSON(rawToChar(x$content))$medianGeneExpression$tissueSiteDetailId,
                                                                        fromJSON(rawToChar(x$content))$medianGeneExpression$median)))

          results_df2 <- lapply(results_df, function(x) dcast(x, V1+V2~V3, value.var = "V4"))
          results_df3 <- lapply(results_df2, function(x) cbind(x))
          GTEx.info <- results_df3[[1]] #initialize the dataframe
          for (i in 2:length(results_df3))
          {
            GTEx.info = bind_rows(GTEx.info, results_df3[[i]])
          }

          GTEx.info[,3:ncol(GTEx.info)] <- apply(GTEx.info[,3:ncol(GTEx.info)], 2, as.numeric)
          colnames(GTEx.info)[c(1,2)] = c("Name", "Description")

          tissue.df <- colMeans(GTEx.info[GTEx.info$Description %in% rownames(mol.unique.df[mol.unique.df[,1] == 1,]),3:ncol(GTEx.info)])
          rownames(mol.unique.df) <- toupper(rownames(mol.unique.df))
          for (i in 2:ncol(mol.unique.df)) {
            tissue.df <- cbind(tissue.df, colMeans(GTEx.info[GTEx.info$Description %in% rownames(mol.unique.df[mol.unique.df[,i] == 1,]),3:ncol(GTEx.info)]))
          }
          colnames(tissue.df) <- colnames(mol.unique.df)
          Object@dfTissue <- tissue.df

          message("[=========================================================]")
          message("[<<<<         ObtainTissueExpression END             >>>>>]")
          message("-----------------------------------------------------------")
          message("[You may want to plot the results using PlotTissueExpression next.]")

          return(Object)
        }
)
