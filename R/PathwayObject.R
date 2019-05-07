#' PathwayObject
#'
#' @slot Data list.
#' @slot PData data.frame.
#' @slot metadata data.frame.
#' @slot plot list.
#' @slot Data.RR data.frame.
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
    Data.RR="data.frame"),
  prototype=list(
    Data=NULL,
    PData=NULL,
    metadata = NULL,
    plot=NULL,
    Data.RR=NULL)
)
