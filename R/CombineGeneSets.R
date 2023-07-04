#' CombineGeneSets
#'
#' Calculate distances between the different experiments.
#'
#' @import clusterProfiler
#' @import limma
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import bigstatsr
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import limma
#' @param Object A PathwayObject.
#' @param combineMethod lets the functions know if the standard RR, a Jaccard index or Cohens Kappa.
#' @param combineMethod.supplied a function which parameter A and parameter B.
#' @param display Either Condensed or Expanded. Expanded seperates the groups
#'                and adds a 1 to every group is this molecular signature is seen here
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
#' man.Great.Object2 <- CombineGeneSets(Object = man.Great.Object1)

#'
#' @export
#'


CombineGeneSets <- function(Object, combineMethod="Standard", combineMethod.supplied, display="Condensed", threads=1)
{
  message("[=========================================================]")
  message("[<<<<            CombineGeneSets START               >>>>>]")
  message("-----------------------------------------------------------")

  #############################################
  ##---------align data stucture ------------##
  #############################################

  if(length(unique(Object@metadata[,"structure"])) > 1)
  {

    message("Warning, data structure not the same, converting all to structure of experiment 1")

    structure.desired <- Object@metadata[1,"structure"]
    Seperator.desired <- as.character(Object@metadata[1,"seperator"])
    for(Meta.i in 2:nrow(Object@metadata))
    {

      gene <- as.character(Object@Data[[Meta.i]][, "Molecules"])
      gene <- strsplit2(x = gene, split = as.character(Object@metadata[Meta.i,"seperator"]))
      for(gene.i in 1:nrow(gene))
      {
        gene.row <-gene[gene.i,]
        gene.row <- gene.row[!gene.row %in% ""]
        gene.df <- bitr(gene.row, fromType = as.character(Object@metadata[Meta.i,"structure"]),
                        toType = c("ENTREZID", "ENSEMBL", "SYMBOL"),
                        OrgDb = as.character(Object@metadata[Meta.i,"organism"]))
        gene.df <- gene.df[,structure.desired]
        Object@Data[[Meta.i]][gene.i, "Molecules"] <- paste(gene.df, collapse = Seperator.desired)

      }
      Object@metadata[Meta.i,"structure"] <- structure.desired
      Object@metadata[Meta.i,"seperator" ] <- Seperator.desired
    }
  }



  #############################################
  ##---------Get unique molecules------------##
  #############################################

  molecules <- vector()
  for(list.i in 1:length(Object@Data))
  {
    x <- strsplit2(as.character(Object@Data[[list.i]]$Molecules), split=as.character(Object@metadata[1,"seperator"]))
    x <- as.vector(x)
    x <- unique(x)
    x <- x[!is.na(x)]
    x <- x[!x ==""]

    molecules <- unique(c(molecules, x))
  }

  #Despite all being symbols, the case might be different, transform everything to the same case
  if(unique(as.character(Object@metadata[,"structure"])) == "SYMBOL")
  {
    message("transforming all genes to upper case, make sure this doesnt change the data")
    message(paste( "raw data has " , length(unique(molecules))," genes", sep=""))
    molecules <- toupper(molecules)
    molecules <- unique(molecules)
    message(paste( "Transformed data has " , length(unique(molecules))," genes", sep=""))
  }


  #############################################
  ##---------combine experiments-------------##
  #############################################

  message("Combining experiments")

  #if its all manually loaded take the data.frame
  if(sum(Object@metadata$loaded == "Manual") == nrow(Object@metadata))
  {
    pathways.i <- as.data.frame(Object@Data)
  }else {
    pathways.i <- do.call(rbind, Object@Data)

  }
  pathways.i2 <- pathways.i[!duplicated(pathways.i$Pathways), ]
  pathways.i$GeneSets <- paste("GeneSet",1:nrow(pathways.i),sep="")

  ################################################
  ##-----------Condense or expand---------------##
  ################################################
  #Condensed or Expanded
  if( display=="Expanded")
  {
    message("preparing expanded display")
    Object@metadata[,"display"] <- rep("Expanded", times = nrow(Object@metadata))

    #If display is expanded then based on the unique molecular signature
    #in which group this signature is in gets marked
    #multiple roups can share the signature
    Pathways.sharedsets <- pathways.i


    for(groups.i in 1:length(Object@PData$Groupnames))
    {
      Pathways.sharedsets$Group <- 0
      colnames(Pathways.sharedsets)[ncol(Pathways.sharedsets)] <- Object@PData$Groupnames[groups.i]
      Pathways.sharedsets[Pathways.sharedsets$Molecules %in% pathways.i[pathways.i$Groups == Object@PData$Groupnames[groups.i], "Molecules"],Object@PData$Groupnames[groups.i]] <- 1
    }
    pathways.i <- Pathways.sharedsets
  }else{
    message("preparing condensed display")

  }

  #####################################
  ##---------Calculate RR------------##
  #####################################
  if(combineMethod == "Standard" | combineMethod == "RR")
  {
    message("calulating RR")

  }else{
    message(paste("calulating",combineMethod))
  }
  pathways.mtx <- matrix(data = 0, nrow = nrow(pathways.i), ncol = length(molecules))
  colnames(pathways.mtx) <- molecules
  rownames(pathways.mtx) <- as.character(pathways.i$Pathways)


  for(list.i in 1:nrow(pathways.i))
  {
    molecules.i <- toupper(as.vector(strsplit2(pathways.i[list.i, "Molecules"], split=Object@metadata[1,"seperator"])))
    pathways.mtx[list.i, molecules.i] <- 1
  }

  ##
  pathways.i$RR_name <- paste(pathways.i[,"Groups"],rownames(pathways.mtx),sep="_")

  ##

  if (threads == 1) {
    RR <- matrix(data = 0, nrow = nrow(pathways.mtx), ncol = nrow(pathways.mtx))
    rownames(RR) <-  pathways.i$RR_name
    colnames(RR) <-  pathways.i$RR_name
  } else {
  RR <- as_FBM(matrix(data = 0, nrow = nrow(pathways.mtx), ncol = nrow(pathways.mtx)))
  }

  ## For RR with unique pathways
  pathways.mtx2 <- matrix(data = 0, nrow = nrow(pathways.i2), ncol = length(molecules))
  colnames(pathways.mtx2) <- molecules
  rownames(pathways.mtx2) <- as.character(pathways.i2$Pathways)

  for(list.i2 in 1:nrow(pathways.i2))
  {
    molecules.i2 <- toupper(as.vector(strsplit2(pathways.i2[list.i2, "Molecules"], split=Object@metadata[1,"seperator"])))
    pathways.mtx2[list.i2, molecules.i2] <- 1
  }

  pathways.i2$RR_name <- rownames(pathways.mtx2)

  if (threads == 1) {
    RR2 <- matrix(data = 0, nrow = nrow(pathways.mtx2), ncol = nrow(pathways.mtx2))
    rownames(RR2) <-  pathways.i2$RR_name
    colnames(RR2) <-  pathways.i2$RR_name
  } else {
    RR2 <- as_FBM(matrix(data = 0, nrow = nrow(pathways.mtx2), ncol = nrow(pathways.mtx2)))
  }

  ############################################
  #---------User supplied function-----------#
  ############################################
  if(combineMethod=="User")
  {
    for(Disease1 in 1:nrow(RR))
    {
      for(Disease2 in 1:ncol(RR))
      {
        RR[Disease1,Disease2] <- combineMethod.supplied(A = pathways.mtx[Disease1,], B = pathways.mtx[Disease2,])
      }
    }
  }

  #############################
  #---------Jaccard-----------#
  #############################
  if(combineMethod=="Jaccard")
  {

    jaccard <- function(A,B)
    {
      M <- sum(as.vector(A) == 1 & as.vector(B) == 1)
      A.c <- sum(as.vector(A) == 1 & as.vector(B) == 0)
      B.c <- sum(as.vector(A) == 0 & as.vector(B) == 1)
      J <- M/(A.c+B.c+M)
      return(J)
    }



    for(Disease1 in 1:nrow(RR))
    {
      for(Disease2 in 1:ncol(RR))
      {
        RR[Disease1,Disease2] <- jaccard(A = pathways.mtx[Disease1,], B = pathways.mtx[Disease2,])
      }
    }
  }

  #############################
  #---------Cohen-----------#
  #############################
  if(combineMethod=="Cohen")
  {

    Cohen <- function(A,B)
    {
      P0 <- ((sum(as.vector(A) == 1 & as.vector(B) == 1)) +
               (sum(as.vector(A) == 0 & as.vector(B) == 0)))/length(A)
      Pe <- ((sum(as.vector(A) == 1) / length(A)) * (sum(as.vector(B) == 1) / length(B))) +
            ((sum(as.vector(A) == 0) / length(A)) * (sum(as.vector(B) == 0) / length(B)))
      K <- ((P0 - Pe)/(1 - Pe))
      return(K)
    }



    for(Disease1 in 1:nrow(RR))
    {
      for(Disease2 in 1:ncol(RR))
      {
        RR[Disease1,Disease2] <- Cohen(A = pathways.mtx[Disease1,], B = pathways.mtx[Disease2,])
      }
    }
  }

  #######################################
  #---------standard function-----------#
  #######################################
  RR.function <- function(x)
  {
    J <- x
    I <- pathways.mtx[Disease1,]
    N <- (ncol(pathways.mtx))#number of molecules

    Pi <- sum(I)#number of molecules in pathway 1
    Pj <- sum(J)#number of molecules in pathway 2

    Cij <- sum(names(J)[J==1] %in% names(I)[I==1])

    if(Cij == 1){Cij <- 0.1}#Adjusted for pathways with 1 molecules
    RRij <- (Cij*N)/((Pi*Pj) - Cij)
    return(RRij)
  }

  RR.function2 <- function(x)
  {
    J <- x
    I <- pathways.mtx2[Disease1,]
    N <- (ncol(pathways.mtx2))#number of molecules

    Pi <- sum(I)#number of molecules in pathway 1
    Pj <- sum(J)#number of molecules in pathway 2

    Cij <- sum(names(J)[J==1] %in% names(I)[I==1])

    if(Cij == 1){Cij <- 0.1}#Adjusted for pathways with 1 molecules
    RRij <- (Cij*N)/((Pi*Pj) - Cij)
    return(RRij)
  }

  if(combineMethod=="Standard")
  {
    if (threads == 1) {
      for(Disease1 in 1:nrow(RR))
      {
        RR[Disease1,] <- apply(X = pathways.mtx, MARGIN = 1, FUN = RR.function )
      }

      for(Disease1 in 1:nrow(RR2))
      {
        RR2[Disease1,] <- apply(X = pathways.mtx2, MARGIN = 1, FUN = RR.function2 )
      }

    } else {#parallelization
        message(paste0("Parallelizing processes, ", threads, " cores have been selected. "), sep="")

        if (threads > detectCores()) {
          message("More cores have been selected than detected.")
          threads <- detectCores()
          message(paste0("The maximum number of cores will be used:", threads, " cores."))
        }

        cl <- makeCluster(threads)
        registerDoParallel(cl)
        foreach(Disease1=1:nrow(RR[])) %dopar% {

          RR[Disease1,] <- apply(X = pathways.mtx, MARGIN = 1, FUN = RR.function )

        }
        stopCluster(cl)

        cl <- makeCluster(threads)
        registerDoParallel(cl)
        foreach(Disease1=1:nrow(RR2[])) %dopar% {

          RR2[Disease1,] <- apply(X = pathways.mtx2, MARGIN = 1, FUN = RR.function2 )

        }
        stopCluster(cl)

        RR <- RR[]
        rownames(RR) <-  pathways.i$RR_name
        colnames(RR) <-  pathways.i$RR_name

        RR2 <- RR2[]
        rownames(RR2) <-  pathways.i2$RR_name
        colnames(RR2) <-  pathways.i2$RR_name
    }
  }
  ######################################
  ##---------Create object------------##
  ######################################

  Object@Data <- list(pathways.i)
  Object@Data.RR <- as.data.frame(RR)
  Object@DataPathways.RR <- as.data.frame(RR2)


  message("-----------------------------------------------------------")
  message("[<<<<<             CombineGeneSets END              >>>>>>]")
  message("[=========================================================]")
  message("[You may want to process ClusterGeneSets next.            ]")


  return(Object)
}
