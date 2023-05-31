#' PlotTreePathway
#'
#' Plots a tree plot with the distances per pathway.
#' @import ggplot2
#' @import tidyverse
#' @import stats
#' @import ggdendro
#' @import cowplot
#' @import clusterProfiler
#' @import ggtree
#' @import patchwork
#' @import RColorBrewer
#'
#' @param Object a pathway object
#' @param clusters A numeric with the number of clusters required
#' @param nodenames boolean to add names to the nodes
#' @param doORA boolean perform ORA analysis
#' @param wordcloud boolean perform boolean to add wordcloud annotations of each cluster
#'
#' @return plot
#'
#' @export
#'
setGeneric(name="PlotTreePathway",
           def=function(Object, clusters = 3, nodenames = TRUE, doORA = TRUE, wordcloud = TRUE)
           {
             standardGeneric("PlotTreePathway")
           }
)


#' PlotTreePathway
#'
#' @param Object a pathway object
#' @param clusters A numeric with the number of clusters required
#' @param nodenames boolean to add names to the nodes
#' @param doORA boolean perform ORA analysis
#' @param wordcloud boolean perform boolean to add wordcloud annotations of each cluster
#'
#' @import ggplot2
#' @import tidyverse
#' @import stats
#' @import ggdendro
#' @import cowplot
#' @import clusterProfiler
#' @import ggtree
#' @import patchwork
#' @import RColorBrewer
#'
#' @return plot
#'
#' @examples

setMethod(f="PlotTreePathway",
          signature = "PathwayObject",
          definition = function(Object, clusters = 3, nodenames = TRUE, doORA = TRUE, wordcloud = TRUE)
          {

          hc <- hclust(as.dist(1 - Object@DataPathways.RR),method = "ward.D")
          # Info of Clusters wanted
          clus <- cutree(hc, clusters)

          keywords <- ORAperCluster(Object, doORA, clusters, clus)

          #For leyend adjust proportionaly to terms length
          keywords_max_length <- max(nchar(as.list(keywords)))
          if (keywords_max_length < 17) {keywords_max_length = 17}

          # Plot dendogram
          g <- split(names(clus), clus)
          p <- ggtree(hc, size=1.2)
          clades <- sapply(g, function(n) MRCA(p, n))
          p <- groupClade(p, clades, group_name='SubTree_ORA') + aes(color=SubTree_ORA)
          mycolors <- c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Set1", n = 6))
          ggtree_plot <- p + scale_color_manual(values=mycolors, breaks=1:clusters, labels = keywords) +
            theme(legend.position='right',legend.justification = c(0,1.5))
          ggtree_plot_noLegend <- ggtree_plot + theme(legend.position = "none")

          log_p.adjust <- -log10(as.numeric(Object@Data[[1]]$Pval))
          dolegend_dot <- TRUE
          if (length(log_p.adjust) == 0) {
            message("P values not detected. Setting all values to 0.05 as default.")
            log_p.adjust <- rep(0.05, length(Object@Data[[1]][,"Pathways"]))
            dolegend_dot <- FALSE
          }

          enrichmentScore <- Object@Data[[1]]$enrichScore
          legendEnrichmentScore <- TRUE
          if (length(enrichmentScore) == 0 | anyNA(Object@Data[[1]]$enrichScore)) {
            message("Enrichment scores not detected.")
            legendEnrichmentScore <- FALSE
          }

          Pathway <- Object@Data[[1]][,"Pathways"]

          if (checkGO(Object)) {
            results <- as.data.frame(unlist(lapply(Pathway, getterm)))
            results$pathway <- Pathway
          } else {
            results <- data.frame(Pathway, Pathway)
          }

          colnames(results) <- c("Term", "Pathway")

          results$Term[which(nchar(results$Term)>50)] <- paste0(substr(results$Term[which(nchar(results$Term)>50)], 1, 30), "...", substr(results$Pathway[which(nchar(results$Term)>50)],1,10), sep="")

          results_unique <- results[!duplicated(results), ]
          rownames(results_unique) <- results_unique$Pathway


          Description <- factor(results$Term,
                                levels = convertNA(results_unique[hc$labels[hc$order], "Term"]))


          Cluster <- factor(Object@Data[[1]][,"Groups"])

          # Prepare for dotplot
          res_plot <- data.frame(Pathway, log_p.adjust, results, Cluster, Description)

          # Plot dotplot
          dotplot <- ggplot(data = res_plot, aes(x=Cluster, y = Description,  size = log_p.adjust)) +
            {if(legendEnrichmentScore==TRUE)aes(color=enrichmentScore)} +
            geom_point() +
            scale_y_discrete(position="right")+
            scale_color_gradient2(low="blue4",mid="white", high="red")+
            theme_cowplot() +
            theme(axis.line = element_blank()) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            {if(nodenames==FALSE)theme(axis.text.y = element_blank())} + #argument show labels
            ylab('') +
            guides(size=guide_legend(title="-log10(p.adj)")) +
            theme(axis.ticks = element_blank(), legend.position = "right", legend.justification = c(0,0))
          dotplot_noLegend <- dotplot + theme(legend.position = "none")

          if (checkGO(Object) == FALSE) {
            message("No GO terms have been detected in the pathways. The semantic enrichment word cloud will not be generated.")
            wordcloud = FALSE
          }

          if (wordcloud == TRUE) {
            plotList <- list()
            for (i in 1:clusters) {
              list <- names(clus[clus == i])
              plot <- wordcloud_generation(list)
              # plot <- plot + annotation_custom(grob = textGrob(label = paste0("Cluster ", i), hjust = 0.5, vjust = -6,
              #                                                  gp = gpar(fontsize = 12, fontface = "bold")))
              plotList[[i]] <- plot + theme(panel.border = element_rect(color = mycolors[i], fill = NA, size = 1))

            }
            word_plot <- plot_grid(plotlist = plotList, ncol = 1)

            } else {
              word_plot <- NULL
          }

          # Legends
          if (doORA == TRUE) {legend_tree <- get_legend(ggtree_plot)} else {legend_tree <- NULL}
          if (dolegend_dot == TRUE) {legend_dot <- get_legend(dotplot)} else {legend_dot <- NULL}

          # Merge Plot
          combine_plot <- plot_grid(ggtree_plot_noLegend, NULL, dotplot_noLegend, nrow= 1, rel_widths= c(0.3,-0.01,2), align = 'h')
          combine_legend <- plot_grid(legend_dot,NULL,legend_tree, ncol=1, rel_heights = c(1,-0.3,1))
          big_plot <- plot_grid(combine_plot, combine_legend, NULL, NULL, word_plot, nrow = 1,
                                rel_widths = c(1, 0.1, keywords_max_length*0.0083, 0.00075, 0.35))

          return(big_plot)
          }
)
