
## User tool: Box plot of metacluster, For a given marker.
## treatmentTable should be a dataframe with two column: "Treatment", "files". Robust specifies either Tukey/lm or Dunn (non adjusted p-values).
## ClustHeat=FALSE for no clustering on heatmap.
BoxPlotMarkerMetaClust = function(TreeMetaCl,Title,treatmentTable,ControlTreatment,BottomMargin,Marker,Robust=TRUE,ClustHeat=TRUE) {
    MarkerIndex = which(TreeMetaCl$fSOMTree$prettyColnames == Marker)
    if (length(MarkerIndex) < 1) {stop("No marker ",Marker," in data")} else {
    BoxPlotMetaClustFull(TreeMetaCl,Title,treatmentTable,ControlTreatment,BottomMargin,yLab="",Norm=FALSE,Marker,Robust,ClustHeat) }
}
