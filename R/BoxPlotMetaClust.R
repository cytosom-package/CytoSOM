## User tool: Box plot of metacluster, either percentage or normalized size is Norm = T.
## TreatmentTable should be a dataframe with two column: "Treatment", "files" (a third one with column "NormalizationFactor" if Norm=T).
## Robust specifies either Tukey/lm or Dunn (non adjusted p-values).
## ClustHeat=FALSE for no clustering on heatmap.
BoxPlotMetaClust = function(TreeMetaCl,Title,treatmentTable,ControlTreatment,BottomMargin,yLab,Norm=FALSE,Robust = TRUE,ClustHeat=TRUE) {
   BoxPlotMetaClustFull(TreeMetaCl,Title,treatmentTable,ControlTreatment,BottomMargin,yLab,Norm,Marker="",Robust,ClustHeat)
}
