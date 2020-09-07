#' ORAGeneSets
#'
#' @import WebGestaltR
#'
#' @importFrom stats dist hclust kmeans
#' @param Object A PathwayObject.
#' @param ORA.returned A numeric for how Go terms should be returned per analysis. Alternative use "All" to get all of them.
#' @param unique.per.cluster A vector with strings of genes within each Gene-Set. Each string is seperated using the seperator supplied.
#' @param check.error If check.error =T then the program will check to make sure the webgestaltR function will run or skip at the cost of preformance
#'
#'
#' @return a data frame with the over represented gene sets for all the clusters.
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
#'ORAGeneSets(Object = IPA.object3,unique.per.cluster = TRUE ,
#'                                              ORA.returned = 1)


ORAGeneSets <- function(Object, ORA.returned = 10, unique.per.cluster=T, check.error=F)
{
  message("[=========================================================]")
  message("[<<<<               ORAGeneSets START                >>>>>]")
  message("-----------------------------------------------------------")

  if(ORA.returned == "All")
  {
    ORA.method = "fdr"
  }else{
    ORA.method = "top"
  }

  clus.mol.ls <- list()
  for(clus.i in 1:max(Object@Data[[1]]$cluster))
  {
    clus.x <- Object@Data[[1]][Object@Data[[1]]$cluster == clus.i,]
    clus.x$Molecules <- as.character(clus.x$Molecules)
    clus.mol.ls[[clus.i]] <- unique(as.vector(strsplit2(clus.x$Molecules, split = Object@metadata$seperator[1])))
  }
  unique.mol <- unique(do.call(what = c, args = clus.mol.ls))

  mol.unique.df <- as.data.frame(matrix(0, nrow = length(unique.mol), ncol = length(clus.mol.ls)))
  rownames(mol.unique.df) <- unique.mol
  colnames(mol.unique.df) <-  paste("Cluster_", 1:length(clus.mol.ls), sep="")
  for(clus.i in 1:max(Object@Data[[1]]$cluster))
  {
    mol.unique.df[clus.mol.ls[[clus.i]],clus.i] <- clus.i
  }

  clus.mol.ls.unique <- list()
  for(clus.i in 1:max(Object@Data[[1]]$cluster))
  {
    idx <- vector()
    for(rows.i in 1:nrow(mol.unique.df))
    {
      if(sum(mol.unique.df[rows.i,] == 0) == (max(Object@Data[[1]]$cluster)-1) & !mol.unique.df[rows.i,clus.i] == 0)
      {
        idx[rows.i] <- T
      }else{
        idx[rows.i] <- F
      }
    }
    clus.mol.ls.unique[[clus.i]] <- rownames(mol.unique.df)[idx]
  }



  Object.list <- list()
  for(clus.i in 1:max(Object@Data[[1]]$cluster))
  {
    
    symbol.type <- Object@metadata[1,"structure"]
    if(symbol.type == "SYMBOL"){symbol.type <- "genesymbol"}

    clus.org <- Object@metadata[1,"organism"]
    if(clus.org == "org.Hs.eg.db"){clus.org <- "hsapiens"}
    if(clus.org == "org.Mm.eg.db"){
      clus.org <- "mmusculus"
    }


    if(unique.per.cluster==F)
    {
      genes.x <- clus.mol.ls[[clus.i]]
    }
    if(unique.per.cluster==T)
    {
      genes.x <- clus.mol.ls.unique[[clus.i]]
    }
    if(length(genes.x) < 3)
    {
      message(paste("Cluster ",clus.i, "has not enough genes",sep=""))
      message("Moving to next cluster")
      next()
    }

    genes.x <- limma::strsplit2(genes.x, split = " ")[,1]
    genes.x <- limma::strsplit2(genes.x, split = "/")[,1]

    genes.x <- genes.x[!genes.x == ""]

    if(check.error == T)
      {
      is.error
function (expr, tell = FALSE, force = FALSE) 
{
    expr_name <- deparse(substitute(expr))
    test <- try(expr, silent = TRUE)
    iserror <- inherits(test, "try-error")
    if (tell) 
        if (iserror) 
            message("Note in is.error: ", test)
    if (force) 
        if (!iserror) 
            stop(expr_name, " is not returning an error.", call. = FALSE)
    iserror
}
       if(is.error(WebGestaltR(enrichMethod = "ORA",
                          organism = clus.org,
                          enrichDatabase = "geneontology_Biological_Process",
                          interestGene =genes.x,
                          interestGeneType  = symbol.type,
                          referenceGeneType = symbol.type,
                          referenceSet = "genome",
                          sigMethod = ORA.method,
                          topThr = ORA.returned,
                          isOutput = F)) == T)
    {
    next()
    }else{
      ORA.clus.i.mol <-  WebGestaltR(enrichMethod = "ORA",
                                        organism = clus.org,
                                        enrichDatabase = "geneontology_Biological_Process",
                                        interestGene =genes.x,
                                        interestGeneType  = symbol.type,
                                        referenceGeneType = symbol.type,
                                        referenceSet = "genome",
                                        sigMethod = ORA.method,
                                        topThr = ORA.returned,
                                        isOutput = F)
                          }
      }else{
    ORA.clus.i.mol <- WebGestaltR(enrichMethod = "ORA",
                                  organism = clus.org,
                                  enrichDatabase = "geneontology_Biological_Process",
                                  interestGene =genes.x,
                                  interestGeneType  = symbol.type,
                                  referenceGeneType = symbol.type,
                                  referenceSet = "genome",
                                  sigMethod = ORA.method,
                                  topThr = ORA.returned,
                                  isOutput = F)
      }
    Object.list[[clus.i]] <- ORA.clus.i.mol
    names(Object.list)[clus.i] <- paste("Cluster_",clus.i,sep="")

  }

  ORA.df <- vector()

  for(clus.i in which(!names(Object.list) == ""))
  {
    x <- Object.list[[clus.i]]
    x$Cluster <- names(Object.list)[clus.i]
    ORA.df <- rbind(ORA.df, x)
  }

  message("-----------------------------------------------------------")
  message("[<<<<<             ORAGeneSets END              >>>>>>]")
  message("[=========================================================]")
  message("[You may want to process HighlightGeneSets next.            ]")
  return(ORA.df)
}
