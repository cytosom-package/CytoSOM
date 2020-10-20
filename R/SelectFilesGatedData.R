#' Select a new data from a list of files
#' @param GatedData data from CytoSOM::DownloadCytoData or CytoSOM::PolygonGatingGatedData
#' @param Files list of file names (with no directory name)
#' @return FlowSOM data object.
#' @export
#'
SelectFilesGatedData <- function(GatedData,Files)
{
  NewRawData = SelectFilesRawData(GatedData$fSOMData,Files)
  NewData = GatedData
  NewData$fSOMData = NewRawData
  return(NewData)
}
  