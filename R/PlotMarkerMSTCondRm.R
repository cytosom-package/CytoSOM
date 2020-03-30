## User tool: marker level represented on metacluster tree, on a subset of samples given by an list of index, removing a given number of smallest metacluster
#' Modified PlotMarkerMST function of FlowSOM, with removal of smallest clusters, over a subset of data
#' @param fSOMObject FlowSOM tree
#' @param markerName name of marker
#' @param condIndex index of data events
#' @param mainTitle title
#' @param nbRm  number of smallest cluster to remove
#' @param globalMinMax min and max value of markers, used for color scale
#' @export
#'
PlotMarkerMSTCondRm <- function(fSOMObject,markerName,condIndex,mainTitle,nbRm,globalMinMax=c()){
    fSOM4Plot=list(
        map=fSOMObject$map,
        prettyColnames=fSOMObject$prettyColnames,
        MST=fSOMObject$MST)
    dataIndex = unlist(sapply(fSOMObject$metaData[condIndex],function(x){x[1]:x[2]}))
    clSizes = sapply(1:length(fSOMObject$map$medianValues[,1]),function(x){length(which(fSOMObject$map$mapping[dataIndex,1] == x))})
    if (nbRm>0) {
        indexKeep =  which(fSOMObject$MST$size > sort(fSOMObject$MST$size)[nbRm])
        indexRemove  = setdiff((1:length(fSOMObject$MST$size)),indexKeep)
        fSOM4Plot$MST$size = (sqrt(clSizes)/max(sqrt(clSizes))*15)[indexKeep]
        fSOM4Plot$map$medianValues = t(sapply(1:length(fSOMObject$map$medianValues[,1]),function(i){
            if(length(intersect(dataIndex,which(fSOMObject$map$mapping[,1] == i))) == 0 )
            {print(paste("Cluster ",i," has size zero for given condition, use global median value"))
                apply(fSOMObject$data[which(fSOMObject$map$mapping[,1] == i),,drop=F],2,function(x){median(x)})}
            else {apply(fSOMObject$data[intersect(dataIndex,which(fSOMObject$map$mapping[,1] == i)),,drop=F],2,function(x){median(x)})}
        }))[indexKeep,]
        ##print(fSOM4Plot$map$medianValues)
        fSOM4Plot$MST$graph=igraph::induced_subgraph(fSOMObject$MST$graph,indexKeep)
        fSOM4Plot$MST$l = fSOMObject$MST$l[indexKeep,]

    }
    else {
        fSOM4Plot$MST$size = sqrt(clSizes)/max(sqrt(clSizes))*15
        fSOM4Plot$map$medianValues = t(sapply(1:length(fSOMObject$map$medianValues[,1]),function(i){
            if(length(intersect(dataIndex,which(fSOMObject$map$mapping[,1] == i))) == 0 )
            {print(paste("Cluster ",i," has size zero for given condition, use global median value"))
                apply(fSOMObject$data[which(fSOMObject$map$mapping[,1] == i),,drop=F],2,function(x){median(x)})}
            else{apply(fSOMObject$data[intersect(dataIndex,which(fSOMObject$map$mapping[,1] == i)),,drop=F],2,function(x){median(x)})}
        }))
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
        FlowSOM::PlotMarker(fSOM4Plot, marker=markerName, view = "MST", main=mainTitle,colorPalette = colorPaletteCond)
     }
     else {FlowSOM::PlotMarker(fSOM4Plot, marker=markerName, view = "MST", main=mainTitle)}
}
