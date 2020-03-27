#' GetSTRINGdbPerGeneSets
#'
#' @import STRINGdb
#'
#' @importFrom stats dist hclust kmeans
#' @param Object A PathwayObject.
#' @param unique.per.cluster A vector with strings of genes within each Gene-Set. Each string is seperated using the seperator supplied.
#' @param plot.input A vector with group names of the different gene set experiments
#'
#'
#' @return a list with the string output and the prep for the plotting function PlotSTRINGdbPerGeneSets
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
#'                                              unique.per.cluster = T ,
#'                                              plot.input = "All")

GetSTRINGdbPerGeneSets <- function(Object, unique.per.cluster=T, plot.input="All")
{

  message("[==================================================================]")
  message("[<<<<             GetSTRINGdbPerGeneSets START                >>>>>]")
  message("--------------------------------------------------------------------")

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
    mol.unique.df[clus.mol.ls[[clus.i]],clus.i] <- 1
  }

  df <- mol.unique.df

  if(Object@metadata$organism[1] == "org.Hs.eg.db")
  {
    species = 9606
  }
  if(Object@metadata$organism[1] == "org.Mm.eg.db")
  {
    species = 10090
  }
  string_db <- STRINGdb$new( version="10", species=species, score_threshold=0, input_directory="" )


  df.plot <- list()
  for(df.i in 1:ncol(df))
  {
    if(unique.per.cluster==F)
    {
      genes.i <- rownames(df[df[,df.i] == 1,])
    }
    if(unique.per.cluster==T)
    {
      genes.i <- rownames(df[df[,df.i] == 1 & rowSums(df) == 1,])
    }
    if(length(genes.i) == 0)
    {
      message(paste("Cluster ",clus.i, "has no unique genes",sep=""))
      message("Moving to next cluster")
      next()
    }

    genes.i <- as.matrix(genes.i)
    colnames(genes.i) <- "gene"
    #print(nrow(genes.i))

    mapped.genes.i <- string_db$map(genes.i, "gene", removeUnmappedRows = TRUE )
    #print(nrow(mapped.genes.i))

    if(plot.input == "All")
    {
      plot.input <- length(mapped.genes.i$STRING_id)
    }

    list.i <- list(map = mapped.genes.i,
                   img = string_db$get_png(mapped.genes.i$STRING_id[1:plot.input], payload_id = NULL, required_score = NULL),
                   input =  string_db$get_summary(mapped.genes.i$STRING_id[1:plot.input]),
                   link =  string_db$get_link(mapped.genes.i$STRING_id[1:plot.input], payload_id = NULL,
                                              required_score = NULL))
    if(unique.per.cluster==T)
    {
      list.i$input <- gsub("proteins", "Unique_GenesPerCluster",list.i$input)
    }

    if(unique.per.cluster==F)
    {
      list.i$input <- gsub("proteins", "All_GenesPerCluster",list.i$input)
    }

    df.plot[[df.i]] <- list.i
  }

  message("-------------------------------------------------------------------")
  message("[<<<<             GetSTRINGdbPerGeneSets ends                >>>>>]")
  message("[=================================================================]")
  message("To view the StringPlots use plotSTRINGdbPerGeneSets")

  return(df.plot)
}


