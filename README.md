# GeneSetCluster

This is a package meant to be used to cluster together Gene-Sets from pathway tools such as Ingenuity Pathway Aanalysis (IPA) (https://www.qiagenbioinformatics.com/products/ingenuity-pathway-analysis/) and GREAT (http://great.stanford.edu/public/html/index.php). Gene-Sets often appear significant when running such tools with different labels displaying different Gene-Sets.

However, it can become difficult to interpret the output sometimes when the data of several gene sets compared. Furthermore the output data has several limitations: 

1) Low ratio: where there are only a few genes enriched. 
2) High similarity: where the different genes sets that appear have the same genes enriched despite the different labels. 
3) Low overlap: where the same gene set labels appear in different experimental settings but different genes are enriched.

So its better to review theses sets of genes together instead of investigating the many different Gene-Sets individually. This package does this by taking the sets of genes of every Gene-Set and calculates the Relative Risk of them appearing together. This is used as a distance between them, the higher it is the greater overlap they have. This distance score is then used to cluster the Gene-Sets together.

The data in the user guide is taken from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111385 Data has been analysed with Deseq2 pipeline where the data is Transcriptome sequencing of WT and conditional-Tgfbr2 knockout microglia and CNS-repopulating monocyte-derived macrophages from C57BL/6 mouse in triplicates. The genes were picked in a comparison of uG vs Mac in both WT and KO with a 1E-06 pvalue cutoff. Genes were analysed in both IPA with the Canonical pathways and functional annotations were exported into a excel file with default settings. Bed files were also generated for the genes, and uploaded in GREAT with a background of the sequencing samples. All data was exported in tsv format.

For how to use the package see https://github.com/TranslationalBioinformaticsUnit/GeneSetCluster/wiki
