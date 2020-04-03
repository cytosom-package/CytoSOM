#' interactive gating from data. Mouse click to build polygon, right click to exit.
#' @param RawData data from FlowSOM::ReadInput
#' @param marker1 x-axis
#' @param marker2 y-axis
#' @param fcsFiles list of file names or list of file indices to be used
#' @param Xlim limits of x-axis, if NULL, limits set automatically
#' @param Ylim lmiits of y-axis, if NULL, limits set automatically
#' @return a polygon with the marker names
#' @export

Interactive2DGating <- function(RawData,marker1,marker2,fcsFiles,Xlim=NULL,Ylim=NULL)
{
  if (typeof(fcsFiles) == "character")
  {
    DataIndices = unlist(lapply(fcsFiles),function(file){
      StartEnd = RawData$metaData[which(gsub(".*/","",names(RawData$metaData)) == file)][[1]];
      (StartEnd[1]:StartEnd[2])
    })
  }
  else
    DataIndices = unlist(lapply(RawData$metaData[fcsFiles],function(StartEnd){(StartEnd[[1]]:StartEnd[[2]])}))
  data2D =
    RawData$data[DataIndices,
                 c(names(which(gsub(" <.*","",RawData$prettyColnames) == marker1)),
                   names(which(gsub(" <.*","",RawData$prettyColnames) == marker2)))]
  tRange = list(c(0,0.15),c(0.015,0.05),c(0.05,0.15))
  myramp = colorRampPalette(c('white', 'red', 'yellow', 'blue', 'green'))

  x11()
  if (is.null(Xlim)) {Xlim=c(min(data2D[,1]),max(data2D[,1]))}
  if (is.null(Ylim)) {Ylim=c(min(data2D[,2]),max(data2D[,2]))}

  suppressWarnings(smoothScatter(data2D[,1],data2D[,2],
                xlab = marker1, ylab = marker2, bandwidth = 0.01,
                nbin=512,nrpoints = 0, colramp = myramp, useRaster = T,xlim=Xlim,ylim = Ylim ))
  CellGate = locator(n=512, type = 'o',col= 'blue')

  return(list(polygon = CellGate, marker1=marker1,marker2=marker2))
}
