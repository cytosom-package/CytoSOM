#' Find the names of metaclusters from a sub-name
#' @param subName sub-part of a meta-cluster name, for meta-cluster number preceeding its name, use "^number_"
#' @param TreeMetaCl FlowSOM tree with meta-clusters, constructed within buildFSOMTree (from the data of FSOMData)
#' @return a vector of names
#' @export

FindMetaClustNames <- function(subName,TreeMetaCl)
  {unique(TreeMetaCl$metaCl)[grep(subName,unique(TreeMetaCl$metaCl),T)]}

