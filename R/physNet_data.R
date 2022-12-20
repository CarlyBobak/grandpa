#' Simulated physician network data
#'
#' Graph simulated from medicare claims wherein
#' physicians are nodes, and edges between nodes indicate overlapping patients.
#' Physician specialty is represented as a vertex attribute and shown using vertex
#' color by default.
#'
#' @docType data
#'
#' @usage data(physNet)
#'
#' @format An object of class \code{"igraph"}; see \code{\link[igraph]{igraph}}.
#'
#' @keywords datasets
#'
#' @references Bobak et al. (2022) arXiv preprint arXiv:2211.15000.
#' (\href{https://arxiv.org/abs/2211.15000}{arXiv})
#'
#' @source \href{https://arxiv.org/abs/2211.15000}{arXiv}
#'
#' @examples
#' data(physNet)
#' plot(physNet,vertex.label=NA)
"physNet"