#' Find the names of metaclusters from a sub-name
#' @param subNames list of sub-part of a meta-cluster name, for meta-cluster number preceeding its name, use "^number_"
#' @param TreeMetaCl FlowSOM tree with meta-clusters, constructed within buildFSOMTree (from the data of FSOMData)
#' @return a vector of names
#' @export

FindMetaClustNames <- function(subNames,TreeMetaCl)
{ if (length(subNames) == 1) {
  unique(TreeMetaCl$metaCl)[grep(gsub("+","\\+",subNames,fixed=T),unique(TreeMetaCl$metaCl),T)]}
  else
  {
  Reduce(intersect,lapply(subNames,function(nm){unique(TreeMetaCl$metaCl)[grep(gsub("+","\\+",nm),unique(TreeMetaCl$metaCl),T)]}))
  }
}
