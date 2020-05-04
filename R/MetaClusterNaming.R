#' Naming of meta-cluster, based on quartile if CutOff is not provided
#' @param TreeMetaCl FlowSOM tree with meta-clusters, contructed for buildFSOMTree
#' @param Markers vector of markers used for naming
#' @param Levels number of level for naming: 2 -> -/+ (1/2-4 quartiles), 3-> -/low/high (1/2-3/4 quartiles), 4-> -/low/int/high (1/2/3/4 quartiles)
#' @param CutOff matrix for optional cut off values for naming: one line for each markers (in the same order), each column for separating levels. If NULL, quartiles are used
#' @param Short is true, low is not represented with level=3, int is not represented if level=4
#' @return data frame with meta-cluster cluster oldName and newName.
#' @export

MetaClusterNaming <- function(TreeMetaCl,Markers,Levels,CutOff=NULL,Short=FALSE)
{
    Levels = as.integer(Levels)
    if ((Levels < 2) | (Levels > 4)) {stop("Not the right amount of levels")}
    MarkerIn = Markers[which(sapply(Markers,function(Marker){length(which(TreeMetaCl$fSOMTree$prettyColnames == Marker))})>0)]
    if (length(setdiff(Markers,MarkerIn))>0) {stop(paste("Don't find markers ",paste(setdiff(Markers,MarkerIn),collapse = " "),sep=""))}
    if (!is.null(CutOff)) {
        if (length(CutOff[,1] != length(Markers))) {stop("not the number of CutOff lines compares to Markers")}
        if (length(CutOff[1,]) != Levels) {steop("number of CutOffs no equal to Levels")}}
    metaClustNewNames=lapply(unique(TreeMetaCl$metaCl),function(metaClust){
        clusterList=which(TreeMetaCl$metaCl == metaClust)
        metaClustIndices=unlist(sapply(clusterList,function(cluster){which(TreeMetaCl$fSOMTree$map$mapping[,1] == cluster)}))
        nameList = sapply(Markers,function(Marker){
            MarkerIndex=which(TreeMetaCl$fSOMTree$prettyColnames == Marker)
            metaClustMedian=median(TreeMetaCl$fSOMTree$data[metaClustIndices,MarkerIndex],na.rm=T)
            if(is.null(CutOff)) {
                vectQuantiles=quantile(TreeMetaCl$fSOMTree$data[,MarkerIndex],na.rm=T)
                if (Levels == 2) {namingMarker = paste(Marker,c("-","+"),sep="")[as.numeric(metaClustMedian > vectQuantiles[2])+1]}
                else if (Levels == 3) {
                    if (Short) {namingMarker = c(paste(Marker,"-",sep=""),"",paste(Marker,"high",sep=""))[
                        as.numeric(metaClustMedian > vectQuantiles[2]) + as.numeric(metaClustMedian > vectQuantiles[4]) +1]}
                    else {namingMarker = c(paste(Marker,"-",sep=""),paste(Marker,"low",sep=""),paste(Marker,"high",sep=""))[
                        as.numeric(metaClustMedian > vectQuantiles[2]) + as.numeric(metaClustMedian > vectQuantiles[4]) +1]}}
                else {if (Short) {namingMarker = c(paste(Marker,"-",sep=""),paste(Marker,"low",sep=""),"",paste(Marker,"high",sep=""))[
                    as.numeric(metaClustMedian > vectQuantiles[2]) +
                        as.numeric(metaClustMedian > vectQuantiles[3]) + as.numeric(metaClustMedian > vectQuantiles[4]) +1]}
                    else {namingMarker = c(paste(Marker,"-",sep=""),paste(Marker,"low",sep=""),paste(Marker,"int",sep=""),paste(Marker,"high",sep=""))[
                        as.numeric(metaClustMedian > vectQuantiles[2]) +
                            as.numeric(metaClustMedian > vectQuantiles[3]) + as.numeric(metaClustMedian > vectQuantiles[4]) +1]}}}
            else {
                cutOffLine = CutOff[which(Markers == Marker),]
                if (Levels == 2) {namingMarker = paste(Marker,c("-","+"),sep="")[as.numeric(metaClustMedian > cutOffLine)+1]}
                else if (Levels == 3) {
                    if (Short) {namingMarker = c(paste(Marker,"-",sep=""),"",paste(Marker,"high",sep=""))[
                        as.numeric(metaClustMedian > cutOffLine[1]) + as.numeric(metaClustMedian > cutOffLine[2]) +1]}
                    else {namingMarker = c(paste(Marker,"-",sep=""),paste(Marker,"low",sep=""),paste(Marker,"high",sep=""))[
                        as.numeric(metaClustMedian > cutOffLine[1]) + as.numeric(metaClustMedian > cutOffLine[2]) +1]}}
                else {if (Short) {namingMarker = c(paste(Marker,"-",sep=""),paste(Marker,"low",sep=""),"",paste(Marker,"high",sep=""))[
                    as.numeric(metaClustMedian > cutOffLine[1]) +
                        as.numeric(metaClustMedian > cutOffLine[2]) + as.numeric(metaClustMedian > cutOffLine[3]) +1]}
                    else {namingMarker = c(paste(Marker,"-",sep=""),paste(Marker,"low",sep=""),paste(Marker,"int",sep=""),paste(Marker,"high",sep=""))[
                        as.numeric(metaClustMedian > cutOffLine[1]) +
                            as.numeric(metaClustMedian > cutOffLine[2]) + as.numeric(metaClustMedian > cutOffLine[3]) +1]}}}
            return(namingMarker)
        })
        nameMeta = paste(nameList,collapse="")
        return(nameMeta)
    })
    metaClustDF=data.frame(oldName=unique(TreeMetaCl$metaCl),newName = unlist(metaClustNewNames))
    return(metaClustDF)
}
