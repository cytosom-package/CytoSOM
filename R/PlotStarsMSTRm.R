#' Plot FlowSOM tree with the capacity to remove the smalles clusters that prevent adequate visualization
#'
#' PlotStarsMSTRm is a modified version of the FlowSOM function PlotStars, that plot the tree build through buildFSOMTree function and allows the removal of a given number of the smallest clusters that could prevent the spatial visualization of the tree.
#' Modified PlotStarsMST function of FlowSOM, with removal of smallest clusters
#' @param fSOMObject FlowSOM tree
#' @param metaClustFactors meta-clusters
#' @param mainTitle title of the plot
#' @param nbRm  number of the smallest clusters to remove
#' @param smallTree true if tree is small and meta-cluster legend is big
#' @param equalSize true if clusters are represented with identical sizes
#' @examples PlotStarsMSTRm(fSOMObject=CytoTree$fSOMTree, metaClustFactors=CytoTree$metaCl,mainTitle="Experiment 1",nbRm=2)
#' @export
#'
PlotStarsMSTRm <- function(fSOMObject,metaClustFactors,mainTitle,nbRm=0,smallTree=F,equalSize=F)
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
       if (equalSize) {
           fSOM4Plot$MST$size = rep(8,length(fSOM4Plot$MST$size))
           }
       PlotStarsBigLeg(fSOM4Plot,backgroundValues = as.factor(metaClustFactors[indexKeep]),
                       main=paste("\n",mainTitle,"\nClusters =",fSOMObject$map$nNodes,", Metaclusters =",length(unique(metaClustFactors)),sep=""),smallTree=smallTree)
    }
    else
    {
        if (equalSize) {fSOM4Plot$MST$size = rep(8,length(fSOM4Plot$MST$size))}
        PlotStarsBigLeg(fSOM4Plot,backgroundValues = as.factor(metaClustFactors),
                    main=paste("\n",mainTitle,"\nClusters =",fSOMObject$map$nNodes,", Metaclusters =",length(unique(metaClustFactors))),smallTree=smallTree)}

}
