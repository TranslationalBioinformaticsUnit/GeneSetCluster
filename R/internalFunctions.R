### Functions --------------------------------------------------------------------
# These functions are not meant to be invoked directly by the user.
# See the PlotPathwayCluster method instead.
obtainDFmetadata <- function(df, mat_cor) {
  ngroups <- unique(df$Groups)
  info <- vector("list", length(ngroups))

  for (i in 1:length(ngroups)) {
    info[[i]] <- df[which(df$Groups==ngroups[i]), "Pathways"]
  }

  df_metadata = as.data.frame(matrix(0, nrow = nrow(mat_cor), ncol = length(ngroups)))
  rownames(df_metadata) = colnames(mat_cor)
  colnames(df_metadata) = ngroups

  for (i in 1:nrow(mat_cor)) {
    for (j in 1:length(ngroups)) {
      if (rownames(mat_cor)[i] %in% info[[j]]) {
        df_metadata[i, j] <- 1
      }
    }
  }

  return(df_metadata)
}


obtainDefCluster <- function(mat_using){

  candidates_clustering <- c()
  definitive_clusters <- list()
  j <- 1
  num_sim = 0.65
  for (i in 1:(ncol(mat_using)-1)) {
    if (i < j-1 | names(mat_using[i,])[j] %in% candidates_clustering) {  ## Avoid repeating same cluster
      next
    }

    if (i == ncol(mat_using)-1) {
      if (mat_using[i,i+1] >= num_sim) {
        candidates_clustering <- c(names(mat_using[i,])[i], names(mat_using[i,])[i+1])
        definitive_clusters <- append(definitive_clusters, list(candidates_clustering))
      }
    } else {

      if (mat_using[i,i+1] >= num_sim) {
        candidates_clustering <- c(names(mat_using[i,])[i])

        for (j in (i+1):(ncol(mat_using))) {

          if (j==ncol(mat_using)) {
            candidates_clustering <- c(candidates_clustering, names(mat_using[i,])[j])
            definitive_clusters <- append(definitive_clusters, list(candidates_clustering))
            break
          }

          if(mat_using[i,j] >= num_sim) {
            candidates_clustering <- c(unique(candidates_clustering), names(mat_using[i,])[j])
            if (j == (ncol(mat_using)-2)) {
              append(definitive_clusters, list(candidates_clustering))
            }

          } else if (mat_using[i+1,j] >= num_sim){ # The lower this value is, the larger group of clusters is obtained
            candidates_clustering <- c(unique(candidates_clustering), names(mat_using[i,])[j])
            if (j == (ncol(mat_using)-2)) {
              append(definitive_clusters, list(candidates_clustering))
            }

          } else {
            definitive_clusters <- append(definitive_clusters, list(candidates_clustering))
            break
          }
        }
      }
    }
  }

  return(definitive_clusters)
}



# This function has been adapted from simplifyEnrichment package (Gu, Z. & Hubschmann, D. 2022).
keywordEnrichment <- function(term_id, tdm, min_bg = 5, min_term = 2) {
  tdm2 = tdm[row_sums(tdm) >= 5, ]

  l = colnames(tdm2) %in% term_id

  n = nrow(tdm2)
  n_term = numeric(n)
  n_bg = numeric(n)
  p = numeric(n)
  for(i in seq_len(n)) {
    if(interactive() && se_opt$verbose) {
      if(i %% 100 == 0 || i == n) {
        message(strrep("\r", 100), appendLF = FALSE)
        message(qq("performing keyword enrichment, @{i}/@{n}..."), appendLF = FALSE)
      }
    }
    v = as.vector(tdm2[i, ])
    s11 = sum(v & l)
    if(s11 < 2) {
      next
    }
    s12 = sum(!v & l)
    s21 = sum(v & !l)
    s22 = sum(!v & !l)

    n_term[i] = s11
    n_bg[i] = s11 + s21

    p[i] = fisher.test(cbind(c(s11, s21), c(s12, s22)), alternative = "greater")$p.value

  }
  if(interactive() && se_opt$verbose) {
    message("")
  }

  df = data.frame(keyword = rownames(tdm2), n_term = n_term, n_bg = n_bg, p = p)
  df = df[df$n_term >= 2, , drop = FALSE]
  df$padj = p.adjust(df$p)
  df[order(df$padj, df$p), , drop = FALSE]
}


#' @import AnnotationDbi
#' @import GO.db
obtainGOdescription <- function(go_id) {
  go_descriptions <- AnnotationDbi::select(GO.db, keys = go_id, columns = "TERM")
  return(go_descriptions)
}


obtainClusterInfo <- function(go_id_list) {

}
