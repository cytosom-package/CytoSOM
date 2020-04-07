#' Apply a list of gates made of polygon (from function Interactive2DGating) to RawData
#' @param GatedData data from CytoSOM::PolygonGatingGatedData
#' @param Polygons list of polygons contructed within Interactive2DGating
#' @param gatingName Name of gating
#' @return FlowSOM data object that includes the polygons gating and the name of the population of interest.
#' @export
#'
PolygonGatingGatedData <- function(GatedData,Polygons,gatingName = "")
{
  NewData <- PolygonGatingRawData(GatedData$fSOMData,Polygons,gatingName)
  if (is.null(GatedData$polygonGating)) {NewPolygonGating = NewData$polygonGating}
  else {NewPolygonGating = c(GatedData$polygonGating,NewData$polygonGating)}
  NewData$polygonGating = NewPolygonGating
  return(NewData)
}
