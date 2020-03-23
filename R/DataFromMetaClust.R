#' Construct a sub-dataset from dataset and a list of metaclusters
#' @param FSOMData FlowSOM data object, created from DownLoadCytoData
#' @param TreeMetaCl tree construct within buildFSOMTree
#' @param MetaClusters vector of metaclusters numbers of names
#' @return a FlowSOM data object
#' @export

DataFromMetaClust <- function(FSOMData,TreeMetaCl,MetaClusters)
{
  newFSOMData=FSOMData[-2]
  Clusters = unlist(lapply(MetaClusters,function(x){which(TreeMetaCl$metaCl == x)}))
  ClusterIndices = unlist(lapply(Clusters,function(x){which(TreeMetaCl$fSOMTree$map$mapping[,1] == x)}))
  newFSOMData$fSOMData$data=FSOMData$fSOMData$data[ClusterIndices,]
  metaDataLengthKept=lapply(FSOMData$fSOMData$metaData,function(x){
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
