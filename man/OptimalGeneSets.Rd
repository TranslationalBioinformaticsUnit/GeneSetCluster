% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OptimalGeneSets.R
\name{OptimalGeneSets}
\alias{OptimalGeneSets}
\title{CombineGeneSets}
\usage{
OptimalGeneSets(object, method, max_cluster, cluster_method, main)
}
\arguments{
\item{object}{A PathwayObject.}

\item{method}{Which method to determing optimal number of clusters. gap, elbow or silhouette.}

\item{max_cluster}{Max number of clusters to test}

\item{cluster_method}{kmeans or hcut. Which clustering method is used}

\item{main}{A string to be used as title in the plot}
}
\value{
a plot
}
\description{
Calculate distances between the different experiments.
}
\examples{

require(GeneSetCluster)
IPA.files <- c(system.file("extdata", "MM10.IPA.KO.uGvsMac.Canonical_pathways.xls",
              package = "GeneSetCluster"),
              system.file("extdata", "MM10.IPA.WT.uGvsMac.Canonical_pathways.xls",
              package = "GeneSetCluster"),
              system.file("extdata", "MM10.IPA.KO.uGvsMac.Functional_annotations.xls",
              package = "GeneSetCluster"),
              system.file("extdata", "MM10.IPA.WT.uGvsMac.Functional_annotations.xls",
              package = "GeneSetCluster"))
canonical.files <- IPA.files[grep("Canonical", IPA.files)]

IPA.object1 <- LoadGeneSets(file_location = canonical.files, #where are  the files
                           groupnames= c("KO", "WT"),
                           P.cutoff = 1.3,
                           Mol.cutoff = 5,
                           Source = "IPA",
                           type = "Canonical_Pathways",
                           structure = "SYMBOL",
                           seperator = ",")
IPA.object2 <- CombineGeneSets(Object = IPA.object1)
OptimalGeneSets(object = IPA.object2, method = "elbow", max_cluster= 24,
               cluster_method = "kmeans", main= "Kmeans for 24 clusters")

}
