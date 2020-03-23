#' Box plot of maker MFI, for the different meta-cluster with statical analysis, output in a PDF file.
#' @param TreeMetaCl FlowSOM tree with meta-clusters, constructed by buildFSOMTree
#' @param Title prefix of the pdf file name
#' @param treatmentTable  data frame containing a column 'files' and a column 'Treatment' (case sensitive)
#' @param Controltreatment name of the control treatment used for statistics
#' @param BottomMargin size of the boxplots bottom margins, in order to be readable
#' @param Marker marker name
#' @param Robust if TRUE, robust statistics is applied (dunn), otherwise linear-model/Tukey
#' @param ClustHeat if TRUE, heatmaps are clustered
#' @return list containing MFI, pairwise comparison p-values, controle vs treatment p-values
#' @export
BoxPlotMarkerMetaClust = function(TreeMetaCl,Title,treatmentTable,
                                  ControlTreatment,BottomMargin,Marker,Robust=TRUE,ClustHeat=TRUE) {
    MarkerIndex = which(TreeMetaCl$fSOMTree$prettyColnames == Marker)
    if (length(MarkerIndex) < 1) {stop("No marker ",Marker," in data")} else {
    BoxPlotMetaClustFull(TreeMetaCl,Title,treatmentTable,ControlTreatment,
                         BottomMargin,yLab="",Norm=FALSE,Marker,Robust,ClustHeat) }
    ## change the name "Size" in "MFI"
}
