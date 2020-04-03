#' Apply a list of gates made of polygon (from function Interactive2DGating) to RawData
#' @param RawData data from FlowSOM::ReadInput
#' @param Polygons list of polygons contructed within Interactive2DGating
#' @param gatingName Name of gating
#' @return FlowSOM data object that includes the polygons gating and the name of the population of interest.
#' @export
#'
PolygonGatingRawData <- function(RawData,Polygons,gatingName = "")
{
  IndexList=lapply(Polygons,function(poly){
    print(poly)
    which(as.logical(sp::point.in.polygon(
      RawData$data[,names(which(gsub(" <.*","",RawData$prettyColnames) == poly$marker1))],
      RawData$data[,names(which(gsub(" <.*","",RawData$prettyColnames) == poly$marker2))],
      poly$polygon$x,poly$polygon$y)))})
  print("Indices constructed")
FullIndex=Reduce(intersect,IndexList)
print("Full Index constructed")
NewData=RawData
NewData$data=RawData$data[FullIndex,]
print("new data constructed")
metaDataLengthKept=lapply(RawData$metaData,function(x){
  length(intersect((x[1]:x[2]),FullIndex))
})
if (length(which(metaDataLengthKept < 1))) {print(paste("Remove files:",names(which(metaDataLengthKept < 1))))}
keepFilesIndices = which(metaDataLengthKept > 0)
LastFilesIndex=cumsum(metaDataLengthKept[keepFilesIndices])
FirstFilesIndex=c(1,LastFilesIndex[-length(LastFilesIndex)]+1)
newMetaData=lapply(1:length(LastFilesIndex),function(x){unname(c(FirstFilesIndex[x],LastFilesIndex[x]))})
names(newMetaData)=names(LastFilesIndex)
print("new metadata constructed")
NewData$metaData=newMetaData
return(list(fSOMData=NewData,polygonGating=Polygons,gatingName = gatingName))
}
