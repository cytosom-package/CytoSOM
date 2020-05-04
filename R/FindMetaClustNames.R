#' Find the names of metaclusters from a sub-name
#' @param subNames list of sub-part of a meta-cluster name, for meta-cluster number preceeding its name
#' @param TreeMetaCl FlowSOM tree with meta-clusters, constructed within buildFSOMTree (from the data of FSOMData)
#' @param start If TRUE, match subnames with the start of meta-cluster names
#' @return a vector of names
#' @export

FindMetaClustNames <- function(subNames,TreeMetaCl,start=F)
{ if (length(subNames) == 1) {unique(TreeMetaCl$metaCl)[grep(subNames,unique(TreeMetaCl$metaCl),fixed=T)]}
  else
  {Reduce(intersect,lapply(subNames,function(nm){
    unique(TreeMetaCl$metaCl)[grep(nm,unique(TreeMetaCl$metaCl),fixed=T)]}))}
  ####
}
