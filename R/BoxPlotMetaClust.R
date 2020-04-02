#' Statistical comparison of metacluster population between experimental condition
#' 
#' BoxPlotMetaClust is a function that allow the statistical analysis, either in percentage or with a normalization, of the impact of experimental condition on each metacluster population.
## TreatmentTable should be a dataframe with two column: "Treatment", "files" (a third one with column "NormalizationFactor" if Norm=T).
## Robust specifies either Tukey/lm or Dunn (non adjusted p-values).
## ClustHeat=FALSE for no clustering on heatmap.
#' Box plot of maker sizes, for the different meta-cluster with statical analysis, output in a PDF file.
#' @param TreeMetaCl FlowSOM tree with meta-clusters, constructed by buildFSOMTree
#' @param Title prefix of the pdf file name
#' @param treatmentTable  data frame containing a column 'files',a column 'Treatment' (case sensitive) and a column "NormalizationFactor" for the optional normalization of the data
#' @param Controltreatment name of the control treatment used for statistics
#' @param BottomMargin size of the boxplots bottom margins, in order to be readable
#' @param ylab name of y-axis of boxplots (eg 'population size')
#' @param Norm if true, sizes are normlized according to column 'NormalizationFactor' of the treatment table. Otherwise, percentage is used
#' @param Robust if TRUE, robust statistics is applied (dunn), otherwise linear-model/Tukey
#' @param ClustHeat if TRUE, heatmaps are clustered
#' @return list containing sizes, pairwise comparison p-values, controle vs treatment p-values
#' @examples BoxplotMetaClust(TreeMetaCl=CytoTree, Title="Experiment1",treatmentTable="ExperimentTable1.csv",ControlTreatment="control",BottomMargin=3,yLab="CD45+",Norm=FALSE,Robust=TRUE,ClustHeat=TRUE)
#' @export
BoxPlotMetaClust = function(TreeMetaCl,Title,treatmentTable,ControlTreatment,BottomMargin,yLab,Norm=FALSE,Robust = TRUE,ClustHeat=TRUE) {
   BoxPlotMetaClustFull(TreeMetaCl,Title,treatmentTable,ControlTreatment,BottomMargin,yLab,Norm,Marker="",Robust,ClustHeat)
}
