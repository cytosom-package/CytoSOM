## User tool: tree representaton of metacluster, given size and marker representation, removing a given number of smallest metacluster
#' Modified PlotStarsMST function of FlowSOM, with removal of smallest clusters
#' @param fSOMObject FlowSOM tree
#' @param metaClustFactors meta-clusters (numbers or names)
#' @param mainTitle title
#' @param nbRm  number of smallest cluster to remove
#' @export
#'
PlotStarsMSTRm <- function(fSOMObject,metaClustFactors,mainTitle,nbRm=0)
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
       PlotStarsBigLeg(fSOM4Plot,backgroundValues = as.factor(metaClustFactors[indexKeep]), main=mainTitle)
    }
    else
        {PlotStarsBigLeg(fSOM4Plot,backgroundValues = as.factor(metaClustFactors), main=mainTitle)}

}
