## ----setup, include=FALSE--------------------------------------------------
require(GeneSetCluster)


## ----Great with Background load--------------------------------------------

Great.files <- c(system.file("extdata", "MM10.GREAT.KO.uGvsMac.bed.tsv", package = "GeneSetCluster"),
                 system.file("extdata", "MM10.GREAT.KO.uGvsMac.bed_BCKGRND.tsv", package = "GeneSetCluster"),
                 system.file("extdata", "MM10.GREAT.WT.uGvsMac.bed.tsv", package = "GeneSetCluster"),
                 system.file("extdata", "MM10.GREAT.WT.uGvsMac.bed_BCKGRND.tsv", package = "GeneSetCluster"))
Great.files.bckgrnd <- Great.files[grepl("BCKGRND", Great.files)]


Great.bckgnrd.Object1 <- LoadGeneSets(file_location = Great.files.bckgrnd, 
                              groupnames= c("KO", "WT"),
                              P.cutoff = 0.05, 
                              Mol.cutoff = 5,
                              Source = "Great",
                              Great.Background = T,#specify the background, as great has a different output if run with or without background
                              type = "Canonical_Pathways",
                              topranks = 20,#Great gives soo much output, recommended is adding a topranks filter for first 20
                              structure = "SYMBOL",
                              Organism = "org.Mm.eg.db",
                              seperator = ",")



## ----Venndiagram, echo=FALSE-----------------------------------------------
par(mfrow=c(2,1))
VennDiagram(n_groups = 2, 
               Group1 = as.character(Great.bckgnrd.Object1@Data$`KO_GO Cellular Component`$Pathways), 
               Group2 =as.character(Great.bckgnrd.Object1@Data$`KO_GO Cellular Component`$Pathways), 
               names_groups = c("KO", "WT"),
               main = "Overlapping GO Cellular Component",legend = F, percentage = F )

VennDiagram(n_groups = 2, 
               Group1 = as.character(Great.bckgnrd.Object1@Data$`KO_GO Biological Process`$Pathways), 
               Group2 =as.character(Great.bckgnrd.Object1@Data$`WT_GO Biological Process`$Pathways), 
               names_groups = c("KO", "WT"),
               main = "Overlapping GO Biological Process",legend = F, percentage = F )

## ----Great with Background meta--------------------------------------------

ShowExperimentdata(Object =Great.bckgnrd.Object1 )
ShowMeta(Object = Great.bckgnrd.Object1 )


## ----Great with Background manage clusters---------------------------------

man.Great.Object1 <- ManageGeneSets(Object = Great.bckgnrd.Object1, keep.type =c("Disease Ontology","GO Biological Process" ), exclude.type="")


ShowExperimentdata(Object =man.Great.Object1 )
ShowMeta(Object =man.Great.Object1 )

## ----Great with Background combine and cluster-----------------------------

man.Great.Object2 <- CombineGeneSets(Object = man.Great.Object1)
man.Great.Object3 <- ClusterGeneSets(Object = man.Great.Object2, 
                                 clusters = 5, 
                                 method = "kmeans")

## ----plotPathways, echo=FALSE----------------------------------------------
PlotGeneSets(Object =man.Great.Object3, fontsize = 3,
            legend = T,
            annotation.mol=F,
            main="Great_Background clustered with Kmeans without scaling \n Disease Ontology and GO Biological Process")

PlotGeneSets(Object =man.Great.Object3, fontsize = 3,
            legend = T,
            annotation.mol=F,
            RR.max = 60,
            main="Great_Background clustered with Kmeans \n Disease Ontology and GO Biological Process")

## ----load IPA--------------------------------------------------------------


IPA.files <- c(system.file("extdata", "MM10.IPA.KO.uGvsMac.Canonical_pathways.xls", package = "GeneSetCluster"),
                 system.file("extdata", "MM10.IPA.WT.uGvsMac.Canonical_pathways.xls", package = "GeneSetCluster"),
                 system.file("extdata", "MM10.IPA.KO.uGvsMac.Functional_annotations.xls", package = "GeneSetCluster"),
                 system.file("extdata", "MM10.IPA.WT.uGvsMac.Functional_annotations.xls", package = "GeneSetCluster"))
canonical.files <- IPA.files[grep("Canonical", IPA.files)]

IPA.object1 <- LoadGeneSets(file_location = canonical.files, #where are  the files
                            groupnames= c("KO", "WT"),#Names of the groups 
                            P.cutoff = 1.3, #minumum cutoff if smaller than 1 it assumes normal pvalue, if larger than 1 it assumes a log10 palue
                            Mol.cutoff = 5,# amount of molecules interested in 
                            Source = "IPA",#How was the data generated
                            type = "Canonical_Pathways",#What is the experiment e.g. canonical pathways, functional anotation
                            structure = "SYMBOL",#structure of the molecules e.g. genenames, ensembl_ID etc 
                            seperator = ",")#How are the genes seperated


## ----Venndiagram IPA Pathways, echo=FALSE----------------------------------
VennDiagram(n_groups = 2, 
               Group1 = as.character(IPA.object1@Data$KO$Pathways), 
               Group2 =as.character(IPA.object1@Data$WT$Pathways), 
               names_groups = c("KO", "WT"),
               main = "Overlapping IPA pathway labels",legend = F, percentage = F )




## ----meta IPA Pathways-----------------------------------------------------
ShowExperimentdata(Object =IPA.object1 )
ShowMeta(Object =IPA.object1 )



## ----combine and cluseter IPA Pathways-------------------------------------

IPA.object2 <- CombineGeneSets(Object = IPA.object1)

IPA.object3 <- ClusterGeneSets(Object = IPA.object2, 
                               clusters = 12, 
                               method = "kmeans")


## ----Highlight IPA Canonical Pathways--------------------------------------

#Highlighting Redox Genes
system.file("data", "Redox.genes.rda", package = "testdat")
IPA.object3.highlight <- HighlightGeneSets(Object = IPA.object3, highligt.genes = Redox.genes, name = "Ros")


## ----PlotPathways with highlight, echo=FALSE-------------------------------

PlotGeneSets(Object =IPA.object3.highlight, fontsize = 3,
            legend = T,
            annotation.mol=F,
            RR.max = 60,
            main="IPA canonical RR with highlights clusters")

## ----User supplied distance calculations-----------------------------------

jaccard <- function(A,B)
{
  #The Jaccard similarity index (sometimes called the Jaccard similarity coefficient) compares members 
  #for two sets to see which members are shared and which are distinct. 
  #It's a measure of similarity for the two sets of data, with a range from 0% to 100%. 
  #The higher the percentage, the more similar the two populations.
  
  M <- sum(as.vector(A) == 1 & as.vector(B) == 1)
  A.c <- sum(as.vector(A) == 1 & as.vector(B) == 0)
  B.c <- sum(as.vector(A) == 0 & as.vector(B) == 1)
  J <- M/(A.c+B.c+M)
  return(J)
}


IPA.Object.J <- CombineGeneSets(Object = IPA.object1, combineMethod = "Jaccard", combineMethod.supplied = jaccard)
IPA.Object.J <- ClusterGeneSets(Object = IPA.Object.J, 
                                                clusters = 4, 
                                                method = "kmeans", 
                                                order = "group")

PlotGeneSets(Object = IPA.Object.J, fontsize =5,
            legend = T,
            annotation.mol=F,
            main="Jaccard distance", RR.max = 50)

## ----load IPA functional annotations---------------------------------------

################################################
#-------IPA on Functional_annotations----------#
################################################


functional.files <- IPA.files[grep("Functional", IPA.files)]

IPA.Functional.object1 <- LoadGeneSets(file_location = functional.files, #where are  the files
                                 groupnames= c("KO", "WT"),#Names of the groups 
                                 P.cutoff = 0.05, #minumum cutoff if smaller than 1 it assumes normal pvalue, if larger than 1 it assumes a log10 palue
                                 Mol.cutoff = 5,# amount of molecules interested in 
                                 Source = "IPA",#How was the data generated
                                 type = "Functional_annotations",#What is the experiment e.g. canonical pathways, Functional_annotations
                                 structure = "SYMBOL",#structure of the molecules e.g. genenames, ensembl_ID etc 
                                 seperator = ",")#How are the genes seperated

## ----Venndiagram IPA functions---------------------------------------------
VennDiagram(n_groups = 2, 
               Group1 = as.character(IPA.Functional.object1@Data$KO$Pathways), 
               Group2 =as.character(IPA.Functional.object1@Data$WT$Pathways), 
               names_groups = c("KO", "WT"),
               main = "Overlapping IPA functions labels",legend = F, percentage = F )




## ----meta IPA functional---------------------------------------------------

ShowExperimentdata(Object =IPA.Functional.object1 )
ShowMeta(Object = IPA.Functional.object1 )



## ----combine and cluseter IPA functional Gene-Sets-------------------------

IPA.Functional.object2 <- CombineGeneSets(Object = IPA.Functional.object1)

IPA.Functional.object3 <- ClusterGeneSets(Object = IPA.Functional.object2, 
                               clusters = 12, 
                               method = "mclust")


## ----plotPathways with functional, echo=FALSE------------------------------

PlotGeneSets(Object = IPA.Functional.object3, 
            fontsize = 3,
            legend = T,
            annotation.mol=F,
            RR.max = 60,
            main="IPA functional RR with mclust clusters")

## ----merge datasets IPA----------------------------------------------------

Ipa.merged.object1 <- MergeObjects(Object.1 =IPA.object1, 
                                    Object.2 = IPA.Functional.object1)



ShowExperimentdata(Object =Ipa.merged.object1 )
ShowMeta(Object = Ipa.merged.object1 )


## ----merge datasets combine and cluster IPA--------------------------------

Ipa.merged.object2 <- CombineGeneSets(Object = Ipa.merged.object1)

Ipa.merged.object3 <- ClusterGeneSets(Object = Ipa.merged.object2, 
                                          clusters = 12, 
                                          method = "kmeans")


## ----highlight dataset merged----------------------------------------------
Ipa.merged.object3.highlight <- HighlightGeneSets(Object = Ipa.merged.object3, highligt.genes = Redox.genes, name = "Ros")


## ----plotPathways with merged, echo=FALSE----------------------------------

PlotGeneSets(Object =Ipa.merged.object3.highlight, fontsize = 3,
            legend = T,
            annotation.mol=F,
            RR.max = 60,
            main="IPA merged RR with highlights clusters")

## ----write merged output to file-------------------------------------------

WriteGeneSets(Object = Ipa.merged.object3.highlight, file_location = "~/Project9/ProjectData_MM10/", name = "Ipa.merged.object3.highlight", write = "Both")

## ----Creating random pathway data------------------------------------------

Test.object <- matrix(data = NA, nrow = 50, ncol = 3)
colnames(Test.object) <- c("Pathways", "Genes", "Group")
Test.object[,"Pathways"] <- paste("Pathway", 1:nrow(Test.object), sep = "_")
Test.object[1:25,"Group"] <- "Group1"
Test.object[26:50,"Group"] <- "Group2"


#Create a random amount of genes per pathway
random.gene <- function()
{
  genenames <- paste("Gene", 1:200, sep = "_")#this gives 200 gene names
  genes <- round(runif(n = runif(n = 1,min =  7,max = 20), min = 1, max = 200), digits = 0)#This gives between 7 and 20 random whole numbers
  genes <- unique(genes)#remove duplicate numbers             
  genes <- genenames[genes]#get random genesnames
  return(genes)
}

for(i in 1:nrow(Test.object))
{
  Test.object[i, "Genes"] <- paste(random.gene(), collapse=",")#Collapse the genes into a string
}

head(Test.object)

## ----creating random pathway Object----------------------------------------

Test.object1 <- ObjectCreator(Pathways = Test.object[,1], 
                              Molecules = Test.object[,2], 
                              Groups = Test.object[,3],
                              Source = "Random",#we randomly generated this data
                              type = "Test",
                              structure = "SYMBOL",
                              sep = ",")#neccesay to seperate the different genes for combinations.
ShowExperimentdata(Object = Test.object1)
ShowMeta(Object = Test.object1)

## ----combining and clustering custom Object--------------------------------
Test.object2 <- CombineGeneSets(Object = Test.object1)
Test.object3 <- ClusterGeneSets(Object = Test.object2, 
                                                clusters = 4, 
                                                method = "kmeans", 
                                                order = "group")


## ----plotting random Object------------------------------------------------
PlotGeneSets(Object = Test.object3, fontsize =5,
            legend = T,
            annotation.mol=F,
            main="Random genes", RR.max = 50)


