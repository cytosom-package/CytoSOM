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
#' @param interQuartile if true, use interquartile (q3-q1) instead of MFI
#' @examples PlotStarsMSTRm(fSOMObject=CytoTree$fSOMTree, metaClustFactors=CytoTree$metaCl,mainTitle="Experiment 1",nbRm=2)
#' @export
#'
PlotStarsMSTRm <- function(fSOMObject,metaClustFactors,mainTitle,nbRm=0,smallTree=F,equalSize=F,interQuartile=F)
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
       if (interQuartile) {
           fSOM4Plot$map$medianValues <-
               t(sapply(indexKeep,
                        function(clust){apply(fSOMObject$data[which(fSOM4Plot$map$mapping[,1] == clust),],2,
                                              function(c){quart = quantile(c);return(quart[4]-quart[2])})}))

       } else
       {fSOM4Plot$map$medianValues=fSOMObject$map$medianValues[indexKeep,]}
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
        if (interQuartile) {
            fSOM4Plot$map$medianValues <-
                t(sapply(1:fSOM4Plot$map$nNodes,
                         function(clust){apply(fSOMObject$data[which(fSOM4Plot$map$mapping[,1] == clust),],2,
                                               function(c){quart = quantile(c);return(quart[4]-quart[2])})}))

        }
        PlotStarsBigLeg(fSOM4Plot,backgroundValues = as.factor(metaClustFactors),
                    main=paste("\n",mainTitle,"\nClusters = ",fSOMObject$map$nNodes,", Metaclusters = ",length(unique(metaClustFactors))),smallTree=smallTree)}

}
