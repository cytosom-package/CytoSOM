
#' build FSOM tree with the metacluster and plot the tree
#' @param fSOMloaded FlowSOM data dowloaded with DowLoadCytoData
#' @param prettyNames list of pretty names needed for building the tree
#' @param clustDim dimension of the 2D cluster grid
#' @param metaClNg number of meta-clusters
#' @param fSOMSeed seed of the random generator used for building the tree
#' @return FlowSOM tree
#' @export
buildFSOMTree <- function(fSOMDloaded,prettyNames,clustDim,metaClNb,fSOMSeed)
{
  set.seed(fSOMSeed)
  if (length(which(names(fSOMDloaded) == "fSOMData")) > 0 )
  {fSOMData = fSOMDloaded$fSOMData} else {fSOMData = fSOMDloaded}
    ## ff<-fSOMDloaded$flJoDataGated$flowSet[[1]]
    fSOMNicePrettyColNames=gsub(" <.*", "", fSOMData$prettyColnames)
    colNamesIndices=unlist(lapply(prettyNames,function(name){which(fSOMNicePrettyColNames == name)}))
    print("Catched col indices:")
    print(colNamesIndices)
    ##channels_of_interest <-  fSOMData$prettyColnames[colNamesIndices]
    fSOM<-FlowSOM::BuildSOM(fSOMData,colsToUse = colNamesIndices,silent = FALSE,xdim=clustDim,ydim=clustDim,rlen=10,init=FALSE,distf=2)
    fSOM<-FlowSOM::BuildMST(fSOM,silent = FALSE,tSNE=FALSE)
    fSOM$prettyColnames =  fSOMNicePrettyColNames
    metacl<-FlowSOM::metaClustering_consensus(fSOM$map$codes,k=metaClNb,seed=fSOMSeed)
    PlotStarsBigLeg(fSOM,backgroundValues = as.factor(metacl))
    return(list(fSOMTree = fSOM,metaCl = metacl))
}
