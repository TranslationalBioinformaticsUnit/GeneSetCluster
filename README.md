# GeneSetCluster

Gene-set analysis (GSA) tools aim to provide biological context by identifying gene-sets associated with a given trait or traits of interest. Gene-set analysis tools are frequently gene-centric, using as input results from studies using RNA-seq or microarrays data (e.g. Ingenuity or GSEA among many others), but it has also been adapted for interval-based analysis derived from data-types such as DNA methylation or ChIP/ATAC-seq. 

Gene-sets are derived from knowledge bases which group genes on the basis of shared biological or functional properties; GSA tools provide, as a result, a list of significant gene-sets. However, while these results are useful for the researcher in the identification of major biological insights, they may be complex to interpret because many gene-sets have largely overlapping gene content or because the result consists of a large number of gene-sets making it complicated to identify the major biological insights. 

To overcome such limitations, we present GeneSetCluster which allows integrating results from one or multiple experiments and/or tools by clustering GSA-derived gene-sets; gene-sets are clustered using a distance metric derived from shared genes between gene-sets. As a result, GeneSetCluster identifies clusters of gene-sets with similar gene-set definitions (i.e. gene content) and, as a result, it allows the researcher to focus on such groups for biological interpretations.

The data in the user guide is taken from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111385 Data has been analysed with Deseq2 pipeline where the data is Transcriptome sequencing of WT and conditional-Tgfbr2 knockout microglia and CNS-repopulating monocyte-derived macrophages from C57BL/6 mouse in triplicates. The genes were picked in a comparison of uG vs Mac in both WT and KO with a 1E-06 pvalue cutoff. Genes were analysed in both IPA with the Canonical pathways and functional annotations were exported into a excel file with default settings. Bed files were also generated for the genes, and uploaded in GREAT with a background of the sequencing samples. All data was exported in tsv format.

# Installation


install.packages("devtools")

require(devtools)

install.github("TranslationalBioinformaticsUnit/GeneSetCluster")


# User guide
For how to use the package see https://github.com/TranslationalBioinformaticsUnit/GeneSetCluster/wiki

![Pipeline](https://github.com/TranslationalBioinformaticsUnit/GeneSetCluster/blob/master/fig/GeneSetCluster_Outline_revision.png)

# Cite this article

Ewing, E., Planell-Picola, N., Jagodic, M. et al. GeneSetCluster: a tool for summarizing and integrating gene-set analysis results. BMC Bioinformatics 21, 443 (2020). https://doi.org/10.1186/s12859-020-03784-z
