#' Select a new data from a list of files
#' @param RawData data from FlowSOM::ReadInput
#' @param Files list of file names (with no directory name)
#' @return FlowSOM data object.
#' @export
#'
SelectFilesRawData <- function(RawData,Files)
{
  IndexList=lapply(Files,function(file){
    (RawData$metaData[which(gsub(".*/","",names(RawData$metaData)) == file)][[1]][1]:
       RawData$metaData[which(gsub(".*/","",names(RawData$metaData)) == file)][[1]][2])
  })
  FullIndex = Reduce(union,IndexList)
  NewData=RawData
  NewData$data=RawData$data[FullIndex,]
  metaDataLengthKept=lapply(RawData$metaData,function(x){
    length(intersect((x[1]:x[2]),FullIndex))
  })
  keepFilesIndices = which(metaDataLengthKept > 0)
  LastFilesIndex=cumsum(metaDataLengthKept[keepFilesIndices])
  FirstFilesIndex=c(1,LastFilesIndex[-length(LastFilesIndex)]+1)
  newMetaData=lapply(1:length(LastFilesIndex),function(x){unname(c(FirstFilesIndex[x],LastFilesIndex[x]))})
  names(newMetaData)=names(LastFilesIndex)
  NewData$metaData=newMetaData
  return(NewData)
}
