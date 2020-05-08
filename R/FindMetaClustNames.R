#' Find each name of metaclusters that contains all sub-names
#' @param subNames list of sub-part of a meta-cluster name
#' @param TreeMetaCl FlowSOM tree with meta-clusters, constructed within buildFSOMTree (from the data of FSOMData)
#' @param start If TRUE, match subnames with the start of meta-cluster names
#' @return a vector of names
#' @export

FindMetaClustNames <- function(subNames,TreeMetaCl,start=F)
{
  if (start) {
    if (length(subNames) == 1) {
      unique(TreeMetaCl$metaCl)[substr(unique(TreeMetaCl$metaCl),start=0,stop=nchar(subNames)) == subNames]
      }
  else
  {
    Reduce(intersect,lapply(subNames,function(nm){
      unique(TreeMetaCl$metaCl)[substr(unique(TreeMetaCl$metaCl),start=0,stop=nchar(nm)) == nm]}))
  }
}
else
{
  if (length(subNames) == 1) {unique(TreeMetaCl$metaCl)[grep(subNames,unique(TreeMetaCl$metaCl),fixed=T)]}
  else
  {Reduce(intersect,lapply(subNames,function(nm){
    unique(TreeMetaCl$metaCl)[grep(nm,unique(TreeMetaCl$metaCl),fixed=T)]}))}
  }
}
