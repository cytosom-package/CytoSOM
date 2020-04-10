#' Construct a sub-dataset from dataset and a list of metaclusters
#' @param FSOMData FlowSOM datobject, created from FlowSOM::ReadInput,
#' DownloadCytoData ($fSOMData element), PolygonGatingRawData ($fSOMData element) or PolygonGatingGatedData ($fSOMData element)
#' @param TreeMetaCl FlowSOM tree with meta-clusters, constructed within buildFSOMTree (from the data of FSOMData)
#' @param MetaClusters vector of metaclusters numbers or names
#' @return FlowSOM data object
#' @export

DataFromMetaClust <- function(FSOMData,TreeMetaCl,MetaClusters)
{
  newFSOMData= list(fSOMData = FSOMData, GateMetaClusters = MetaClusters)
  Clusters = unlist(lapply(MetaClusters,function(x){which(TreeMetaCl$metaCl == x)}))
  ClusterIndices = unlist(lapply(Clusters,function(x){which(TreeMetaCl$fSOMTree$map$mapping[,1] == x)}))
  newFSOMData$fSOMData$data=FSOMData$data[ClusterIndices,]
  metaDataLengthKept=lapply(FSOMData$metaData,function(x){
    length(intersect((x[1]:x[2]),ClusterIndices))
  })
  if (length(which(metaDataLengthKept < 1))) {print(paste("Remove files:",names(which(metaDataLengthKept < 1))))}
  keepFilesIndices = which(metaDataLengthKept > 0)
  LastFilesIndex=cumsum(metaDataLengthKept[keepFilesIndices])
  FirstFilesIndex=c(1,LastFilesIndex[-length(LastFilesIndex)]+1)
  newMetaData=lapply(1:length(LastFilesIndex),function(x){unname(c(FirstFilesIndex[x],LastFilesIndex[x]))})
  names(newMetaData)=names(LastFilesIndex)
  newFSOMData$fSOMData$metaData=newMetaData
  return(newFSOMData)
}
