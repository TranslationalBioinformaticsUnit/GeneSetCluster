setClassUnion("DForNULL", c("data.frame", "NULL", "matrix"))

#' PathwayObject
#'
#' @slot Data list.
#' @slot PData data.frame.
#' @slot metadata data.frame.
#' @slot plot list.
#' @slot Data.RR data.frame.
#' @slot DataPathways.RR data.frame.
#' @slot dfTissue data.frame.
#'
#' @return a PathwayObject
#' @export
#'
#'
#'

PathwayObject <- setClass(
  Class="PathwayObject",
  slots=list(
    Data="list",
    PData="data.frame",
    metadata="data.frame",
    plot="list",
    Data.RR="data.frame",
    DataPathways.RR="DForNULL",
    dfTissue="DForNULL"),
  prototype=list(
    Data=NULL,
    PData=NULL,
    metadata = NULL,
    plot=NULL,
    Data.RR=NULL,
    DataPathways.RR=NULL,
    dfTissue=NULL)
)
