#' ObjectCreator
#'
#' Manually create a PathwayObject
#'
#' @param Pathways A vector with Gene-Set labels.
#' @param Molecules A vector with strings of genes within each Gene-Set. Each string is seperated using the seperator supplied.
#' @param Groups A vector with group names of the different gene set experiments
#' @param Source Tool used to generate gene sets..
#' @param type For IPA data if Canonical pathways or Functional Anotations were supplied.
#' @param structure The structure of the genes. is it SYMBOLS, ENSEMBL, NCBI etc. Used for converting when there is mutiple structure in the object.
#' @param Organism the package name for the human or mouse data, used for converting the gene structure. name of the package, currently org.Hs.eg.db and org.Mm.eg.db supported.
#' @param seperator A character used within in the string to seperate genes
#'
#' @return a pathwayobject
#'
#' @export
#' @examples
#' Test.object <- matrix(data = NA, nrow = 50, ncol = 3)
#'colnames(Test.object) <- c("Pathways", "Genes", "Group")
#'Test.object[,"Pathways"] <- paste("Pathway", 1:nrow(Test.object),
#'                                  sep = "_")
#'Test.object[1:25,"Group"] <- "Group1"
#'Test.object[26:50,"Group"] <- "Group2"
#'
#'
#'#Create a random amount of genes per pathway
#'random.gene <- function()
#'{
#'  genenames <- paste("Gene", 1:200, sep = "_")#this gives 200 gene names
#'  genes <- round(runif(n = runif(n = 1,min =  7,max = 20),
#'                       min = 1, max = 200), digits = 0)
#'  #This gives between 7 and 20 random whole numbers
#'  genes <- unique(genes)#remove duplicate numbers
#'  genes <- genenames[genes]#get random genesnames
#'  return(genes)
#'}
#'
#'for(i in 1:nrow(Test.object))
#'{
#'  Test.object[i, "Genes"] <- paste(random.gene(), collapse=",")
#'}
#'Test.object1 <- ObjectCreator(Pathways = Test.object[,1],
#'Molecules = Test.object[,2],
#'Groups = Test.object[,3],
#'Source = "Random",#we randomly generated this data
#'type = "Test",
#'structure = "SYMBOL",
#'sep = ",")#neccesay to seperate the different genes for combinations.
#'
#'
#'

#'
ObjectCreator <- function(Pathways, Molecules, Groups, Source, Type, structure, sep, organism = NA)
{
  message("[=========================================================]")
  message("[<<<<            ObjectCreator START                 >>>>>]")
  message("-----------------------------------------------------------")
  ###########################
  #--------Data-----------#
  ###########################

  Data <- matrix(NA, nrow=length(Pathways), ncol = 8)
  colnames(Data) <- c("Pathways", "Molecules","Type", "Groups", "Pval", "Ratio", "MoleculesCount", "GeneSet")
  Data <- as.data.frame(Data)

  Data$Pathways <- Pathways
  Data$Molecules <- Molecules
  Data$Groups <- Groups
  Data$Type <- rep(Type, times = nrow(Data))

  for(can.i in 1:nrow(Data))
  {
    Data[can.i,"MoleculesCount"]<- length(as.vector(strsplit2(as.character(Data[can.i,"Molecules"]), split=sep)))
  }

  ###########################
  #--------Pdata-----------#
  ###########################
  Pdata <- matrix(NA,nrow=length(unique(Groups)), ncol = 3)
  colnames(Pdata) <- c("Groupnames", "Length", "file_location")
  rownames(Pdata) <- paste("Experiment",1:length(unique(Groups)), sep="_" )
  Pdata <- as.data.frame(Pdata)

  Pdata$Groupnames <- unique(Groups)
  Pdata[,"Length"] <- table(Groups)

  ###########################
  #--------metadata-----------#
  ###########################
  metadata <- matrix(NA, nrow = length(unique(Groups)), ncol = 12)
  colnames(metadata) <- c("source", "type", "structure", "organism", "seperator", "Data",
                          "cluster.method", "highlight", "order.group", "loaded", "display", "mol.signature")
  rownames(metadata) <- paste("Experiment",1:length(unique(Groups)), sep="_" )
  metadata <- as.data.frame(metadata)

  metadata$structure <- structure
  metadata$organism <- organism

  metadata$source <- Source
  metadata$type <- Type
  metadata$seperator <- rep(sep, times = nrow(metadata))
  metadata$Data <- rep("DF", times = nrow(metadata))
  metadata$loaded <- rep("Manual", times = nrow(metadata))
  metadata[,"display"] <- rep("Condensed", times = nrow(metadata))
  metadata[,"mol.signature"] <- rep("All", times = nrow(metadata))

  ###########################
  #--------Return-----------#
  ###########################
  Object <-  PathwayObject(Data = list(Data),
                           PData = Pdata,
                           metadata = metadata,
                           Data.RR = data.frame(),
                           plot = list(aka2 = data.frame(),
                                       aka3 = data.frame()))


  message("-----------------------------------------------------------")
  message("[<<<<<               ObjectCreator END              >>>>>>]")
  message("[=========================================================]")
  message("[You may want to process CombineGeneSets next.            ]")
  message("[or merge objects using MergeObjects                     ]")
  message("[or select certain types from objects using ManageGeneSets]")

  return(Object)
}
