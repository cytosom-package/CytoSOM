## User tool: marker level represented on metacluster tree, removing a given number of smallest metacluster
#' Modified PlotMarkerMST function of FlowSOM, with removal of smallest clusters
#' @param fSOMObject FlowSOM tree
#' @param markerName name of marker
#' @param mainTitle title
#' @param nbRm  number of smallest cluster to remove
#' @param globalMinMax min and max value of markers, used for color scale
#' @export
#'
PlotMarkerMSTRm <- function(fSOMObject,markerName,mainTitle,nbRm=0,globalMinMax=c())
{
   fSOM4Plot=list(
        map=fSOMObject$map,
        prettyColnames=fSOMObject$prettyColnames,
       MST=fSOMObject$MST)
    if (nbRm>0)
    {
       indexKeep =  which(fSOMObject$MST$size > sort(fSOMObject$MST$size)[nbRm])
       indexRemove  = setdiff((1:length(fSOMObject$MST$size)),indexKeep)
       fSOM4Plot$MST$size=fSOMObject$MST$size[indexKeep]
       fSOM4Plot$map$medianValues=fSOMObject$map$medianValues[indexKeep,]
       fSOM4Plot$MST$graph=igraph::induced_subgraph(fSOMObject$MST$graph,indexKeep)
       fSOM4Plot$MST$l = fSOMObject$MST$l[indexKeep,]
    }

 if (length(globalMinMax) > 0)
     {
         maxGlobal = globalMinMax[2]
         minGlobal = globalMinMax[1]
         maxCond=max(fSOM4Plot$map$medianValues[,markerName])
         minCond=min(fSOM4Plot$map$medianValues[,markerName])
         colorPalette1000=grDevices::colorRampPalette(c("#00007F", "blue","#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(1000)
        minIndex=(minCond-minGlobal)/(maxGlobal-minGlobal)*1000
        maxIndex=(maxCond-minGlobal)/(maxGlobal-minGlobal)*1000
        colorIndexRound=round(minIndex+(0:8)*(maxIndex-minIndex)/8)
        colorPaletteCond=grDevices::colorRampPalette(colorPalette1000[colorIndexRound])

         FlowSOM::PlotMarker(fSOM4Plot,marker=markerName, view = "MST",main=mainTitle,colorPalette = colorPaletteCond)
     }
   else {FlowSOM::PlotMarker(fSOM4Plot,marker=markerName, view = "MST",main=mainTitle)
   }
}
