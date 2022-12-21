


ORA2Cluster.List.fc <- function(Data,
                                structure = "genesymbol",
                                enrichMethod = "ORA", 
                                organism = "org.Hs.eg.db",
                                enrichDatabase = "geneontology_Biological_Process",
                                Mol.cutoff = 5)
{
  message("[=========================================================]")
  message("[<<<<               ORA2Cluster.List.fc              >>>>>]")
  message("-----------------------------------------------------------")
  
  require(WebGestaltR)
  if(is.list(Data) ==F)
  {
    message("Data is not a list")
    stop("Exiting function")
  }
  
  
  ########################################################################3
  message("Identifying Genes")
  Genes.id.ls <- list()
  for(i in 1:length(Data))
  {
    Genes.id.ls[[i]] <- colnames(Data[[i]]) %in% c("Gene", "gene", "genes", "Genes")
    if(sum(Genes.id.ls[[i]]) == 0)
    {
      message("Genes not identified")
      stop("Exiting function")
    }
  }
  
  #########################################################################
  if(structure == "SYMBOL"){structure <- "genesymbol"}
  
  if(organism == "org.Hs.eg.db"){organism2 <- "hsapiens"}
  if(organism == "org.Mm.eg.db" | organism == "mm10" | organism == "MMusculus"){
    organism2 <- "mmusculus"
  }
  message("Starting ORA")
  
  ORA.ls <- list()
  for(i in 1:length(Data))
  {
    message(paste0("ORA on ", names(Data)[i]))
    
    genes.x <- Data[[i]][,Genes.id.ls[[i]]]
    genes.x <- unique(genes.x)
    genes.x <- genes.x[!genes.x == ""]
    message(paste0("Identified ", length(genes.x), " genes"))
    temp.df <-  WebGestaltR(enrichMethod = enrichMethod,
                            organism = organism2,
                            enrichDatabase = enrichDatabase,
                            interestGene =as.character(genes.x),
                            interestGeneType  = structure,
                            referenceGeneType = structure,
                            referenceSet = "genome",
                            sigMethod = "fdr",
                            fdrMethod = "BH",
                            fdrThr = 0.05,
                            isOutput = F)
    temp.df$Group <- names(Data)[i]
    ORA.ls[[i]] <- temp.df
  }
  if(length(ORA.ls) == 0)
  {
    message("No significant pathways identified")
    stop("Exiting function")
  }
  ####################################################
  ####################################################
  #------------Making an GSC object--------------#####
  ####################################################
  ####################################################
  
  message("Making an GSC object")
  
  ORA.df <- do.call(rbind, ORA.ls)
  
  x <- strsplit2(ORA.df$userId, split=";")
  ORA.df$nMolecules <- ncol(x) - rowSums(x == "")
  ORA.df <- ORA.df[ORA.df$nMolecules > Mol.cutoff,]
  
  if(structure == "genesymbol")
  {
    structure <- "SYMBOL"
  }
  
  Object <- ObjectCreator(Pathways= paste0(ORA.df$geneSet, "_", ORA.df$description), 
                          Molecules = ORA.df$userId, 
                          Groups = ORA.df$Group, 
                          Source = enrichMethod, 
                          Type = enrichDatabase, 
                          structure = structure, 
                          sep = ";", 
                          organism = organism)
  
  Object@Data[[1]]$Pval <- ORA.df$FDR
  Object@Data[[1]]$Ratio <- ORA.df$enrichmentRatio
  
  message("-----------------------------------------------------------")
  message("[<<<<<               ORA2Cluster.List.fc END         >>>>>>]")
  message("[=========================================================]")
  message("[You may want to process CombineGeneSets next.            ]")
  message("[or merge objects using MergeObjects                     ]")
  message("[or select certain types from objects using ManageGeneSets]")
  return(Object)
}


#Data <- list(L1 = PBMC.RTX.HCvsBS_BSvsBS.Age.Sex.Duration$Contrasts$HCvsB1.BS$limma[PBMC.RTX.HCvsBS_BSvsBS.Age.Sex.Duration$Contrasts$HCvsB1.BS$limma$P.Value < 0.001,],
#             L2 = PBMC.RTX.HCvsBS_BSvsBS.Age.Sex.Duration$Contrasts$HCvsB2.BS$limma[PBMC.RTX.HCvsBS_BSvsBS.Age.Sex.Duration$Contrasts$HCvsB2.BS$limma$P.Value < 0.001,])
#

#Data$L1$Gene <- as.character(probe.features.merged[rownames(Data$L1),"gene"])
3Data$L2$Gene <- as.character(probe.features.merged[rownames(Data$L2),"gene"])
#
#Test.Object <-ORA2Cluster.List.fc (Data,
#                     structure = "genesymbol",
#                     enrichMethod = "ORA", 
#                     organism = "org.Hs.eg.db",
#                     enrichDatabase = "geneontology_Biological_Process",
#                     Mol.cutoff = 5)
