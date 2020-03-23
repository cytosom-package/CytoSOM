#' Box plot of maker MFI, for the different meta-cluster, in a PDF file.
#' @param TreeMetaCl a FlowSOM tree
#' @param Title Prefix of the pdf file name
#' @param treatmentTable  data frame containing a column 'files' and a column 'Treatment' (case sensitive)
#' @param Controltreatment the name of the control treatment used for statistics
#' @param BottomMargin size of the boxplots bottom margins, in order to be readable
#' @param Marker marker name
#' @param Robust if TRUE, robust statistics is applied (dunn), otherwise linear-model/Tukey
#' @param ClustHeat if TRUE, heatmaps are clustered
#' @return a list containing MFI, pairwise comparison p-values, controle vs treatment p-values
#' @export
BoxPlotMarkerMetaClust = function(TreeMetaCl,Title,treatmentTable,
                                  ControlTreatment,BottomMargin,Marker,Robust=TRUE,ClustHeat=TRUE) {
    MarkerIndex = which(TreeMetaCl$fSOMTree$prettyColnames == Marker)
    if (length(MarkerIndex) < 1) {stop("No marker ",Marker," in data")} else {
    BoxPlotMetaClustFull(TreeMetaCl,Title,treatmentTable,ControlTreatment,
                         BottomMargin,yLab="",Norm=FALSE,Marker,Robust,ClustHeat) }
    ## change the name "Size" in "MFI"
}
