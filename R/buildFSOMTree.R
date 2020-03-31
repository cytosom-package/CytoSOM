#' Build and plot FlowSOM tree 
#' 
#' BuildFSOMTree is a function that plot a tree based on the FlowSOM clusterization and metaclusterization algorithm. From data of interest contained in the DownLoadCytoData object, buildFSOMTree applies FlowSOM algorithm to cluster and metacluster data thanks to the markers of interest.   
#' @param fSOMDloaded FlowSOM data dowloaded with DowLoadCytoData
#' @param prettyNames list of markers names to include to build the tree
#' @param clustDim dimension of the 2D cluster grid
#' @param metaClNb number of metaclusters
#' @param fSOMSeed seed of the random generator used for building the tree
#' @return FlowSOM object that includes the correspondance between the clusters and the metaclusters and the name of the population of interest.
#' @examples CytoTree=buildFSOMTree(fSOMloaded=CytoData,prettyNames=c("CD4","CD8","FoxP3"),ClustDim=4,metaClNb=6,fSOMSeed=1)
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
    if (is.null(fSOMDloaded$gatingName)){gName = ""} else {gName = fSOMDloaded$gatingName}
    return(list(fSOMTree = fSOM,metaCl = metacl,gatingName=gName))
}
