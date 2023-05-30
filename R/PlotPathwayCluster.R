#' PlotPathwayCluster
#'
#' Plots a correlation matrix showing the correlation between the pathways and a functional annotation of the group of pathways based on GO terms descriptions.
#' @import simplifyEnrichment
#' @import slam
#' @import GetoptLong
#' @import ComplexHeatmap
#' @import grid
#'
#' @param Object a pathway object
#' @param nPathway minimum number of pathways per group
#' @param doORA boolean perform ORA analysis
#'
#' @return plot
#'
#' @export
#'
setGeneric(name="PlotPathwayCluster",
           def=function(Object, nPathway = 6, doORA = TRUE, wordcloud = TRUE)
           {
             standardGeneric("PlotPathwayCluster")
           }
)

#' PlotPathwayCluster
#'
#' @param Object a pathway object
#' @param nPathway minimum number of pathways per group
#' @param doORA boolean perform ORA analysis
#'
#' @return plot
#'
#' @examples

setMethod(f="PlotPathwayCluster",
          signature = "PathwayObject",
          definition = function(Object, nPathway = 6, doORA= TRUE, wordcloud = TRUE)
{
  df <- Object@DataPathways.RR
  maxgo <- max(df)
  ordergo <- hclust(1 - as.dist(df))$order
  mat_cor <- cor(df[ordergo, ordergo])

  res <- obtainDefCluster(mat_cor)

  if (length(res[which(lapply(res, function(x) length(x))>=nPathway)]) == 0) {
    message(paste0("There is no groups with at least ", nPathway, " pathways. Using the minimum value, 2 pathways per group."))
    go_id_list <- res[which(lapply(res, function(x) length(x))>=2)]
  } else {
    go_id_list <- res[which(lapply(res, function(x) length(x))>=nPathway)]
  }

  #ora
  clusters <- length(go_id_list)
  clus <- vector()
  for (i in 1:clusters) {
    paths <- c(rep(i, length(go_id_list[[i]])))
    names(paths) <- go_id_list[[i]]
    clus <- c(clus, paths)
  }

  keywords_ora <- ORAperCluster(Object, doORA, clusters, clus)

  if (checkGO(Object) == FALSE) {
    message("No GO terms have been detected in the pathways. The semantic enrichment word cloud will not be generated.")
    wordcloud = FALSE
  } else {
    #adapted from anno_word_cloud_from_GO function
    env_tdm_GO <- readRDS(system.file("extdata", "tdm_GO.rds", package = "simplifyEnrichment"))
    names(go_id_list) <- as.character(1:length(go_id_list))

    # keyword enrichment
    message(paste0("Performing keyword enrichment for "), length(go_id_list), " group(s) of pathways.")
    term <- lapply(go_id_list, function(x, min_stat=5) {
      df <- keywordEnrichment(x, env_tdm_GO)
      df <- df[df$p <= min_stat, , drop = FALSE]
      data.frame(df[, 1], -log10(df$p))
    })
  }

  align_to <- go_id_list

  for (i in 1:length(align_to)){
    for (j in 1:length(align_to[[i]]))
      align_to[[i]][j] <- as.numeric(which(colnames(mat_cor) == align_to[[i]][j]))
  }

  align_to <- lapply(align_to, as.numeric)


  if (wordcloud == TRUE) {
    annot_label <- rowAnnotation(keywords = anno_word_cloud(align_to, term),
                                 annotation_name_align=T)
  } else {
    annot_label <- NULL
  }

  rowNumbers <- rowAnnotation(foo = anno_mark(at = unlist(lapply(align_to, function(x) round(mean(x)))),
                                    labels = rep(paste0("Group ", 1:length(align_to))),
                                    link_width = unit(0,"mm")))

  df_metadata <- obtainDFmetadata(Object@Data[[1]], mat_cor)

  annot_top <- HeatmapAnnotation(df=df_metadata, show_legend = F, annotation_name_gp = gpar(fontsize = 10))
  ht_opt$message <- FALSE
  plot <- Heatmap(mat_cor, cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F ,
                  name="Similarity", top_annotation = annot_top, right_annotation = annot_label) + rowNumbers

  if (doORA == TRUE) {
    legend <- Legend(labels = rep(paste0(1:length(keywords_ora), "- ", keywords_ora)),
                     title = "\nORA", legend_gp = gpar(fill = 1:length(keywords_ora)),
                     nr=8,  title_position = "leftcenter-rot")

    mergedplot <- draw(plot, annotation_legend_list=legend, annotation_legend_side = "bottom")
  } else {
    mergedplot <- plot
  }

  return(mergedplot)
}
)

