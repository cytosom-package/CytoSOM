CytoSOM: a package for easy use of FlowSOM
==========================================
Introduction
==============
CytoSOM helps to apply [FlowSOM](https://github.com/SofieVG/FlowSOM) to flow cytometry data, and perform statistical analysis.

## Installation

```R
devtools::install_github(repo ="gautierstoll/CytoSOM")
```

## HowTo

`CytoSOM` is applied when data has been gated within FlowJo software. A gating strategy without FlowJo is under construction. 

`CytoSOM` needs to be installed. In a folder (eg `MyData`), put your FlowJo sesssion file (eg `file.wsp`) and a sub-folder (eg `FCSdata`) that contains the `.fcs` files (eg `Tube`...`.fcs`)

1. Launch `RStudio`

2. Set you working direcory as the folder where your FlowJo session and you data are, using the full folder name:
```R
setwd("full_name_before_myData/MyData")
```

3. Create an R object that collect your data (eg `CytData`), indicating the name of the gating used in FlowJo (eg `CD45`):
```R
CytoData <- CytoSOM::DownLoadCytoData(dirFCS="FCSdata","CD45",fcsPattern = "Tube",compensate=FALSE )
```

4. Build a clustering tree (eg `CytoTree`), indicating the list of makers used for cluster (eg `c("CD4","CD8","CD11b","FOXP3")`), the size of the cluster grid (eg `7`), the number of meta-clusters (eg `10`), and the seed of the random generator (eg `0`):
```R
CytoTree <- CytoSOM::buildFSOMTree(CytoData,c("CD3","CD4","CD8","CD11b","FOXP3", "CD19"),7,25,1)
```

5. Plot the tree
```R
CytoSOM::PlotStarsMSTRm(CytoTree$fSOMTree,CytoTree$metaCl,"Title Name",0)
```

If the smallest clusters make the image diffcult ot interpret, you can plot the tree by removing the smallest clsuters (eg `3`)
```R
CytoSOM::PlotStarsMSTRm(CytoTree$fSOMTree,CytoTree$metaCl,"Title Name",3)
```

If the tree look staitfactory, you can do the statsistical analysis. Otherwise, try to rebuild the tree with other paramters (cluster grid size and/or meta-cluster number)

6. Rename the meta-cluster names, using a set of markers (eg `c("CD4","CD8","CD11b","FOXP3", "CD19")`)
```R
CytoTree <- CytoSOM::TreeMetaRenaming(CytoTree,c("CD4","CD8","CD11b","FOXP3", "CD19"),"shortRobustName")
```

7. Download an annotation table in `.csv` format (eg `AnnotationTable.csv`), that contains a column 'files', a column 'Treatment' and an optinal column 'Normalization factor', indiacting the separator (eg `;`):

```R
tableTreatmentFCS <- read.csv("AnnotationTable.csv")
```

8. Plot trees in `.pdf` file, for different population sizes and different markers (eg `c("MHCII","PD1","PDL1","PDL2")`):
```R
CytoSOM::plotTreeSet(CytoTree ,c("MHCII","PD1","PDL1","PDL2"),"TitleTrees",rmClNb=0,tableTreatmentFCS,globalScale=T)
```
This produces two pdf in the working directory: `TitleTrees_ClusterTree.pdf` and `TitleTrees_MarkerTree.pdf`.

9. Perform statistical analysis of meta-cluster sizes:
```R
StatAnalysisSizes <-CytoSOM:: BoxPlotMetaClust (CytoTree,Title = "MyTitle",tableTreatmentFCS,ControlTreatmen="PBS",BottomMargin=3,yLab="CD45",Norm=FALSE,Robust = TRUE,ClustHeat=TRUE)
```
If `Norm` is set to `TRUE`, the column 'NormalizationFactor' is used to normalize the meta-cluster sizes. Otherwise, the analisis is perfomed on relative size (precentage). A file `MyTitle_BoxPlotPercentMetaClust.pdf` is produced. Two files containing p-values are also produced: `MyTitle_LmPvalPercentMetacl.csv` and `MyTitle_PairwisePvalPercentMetacl.csv`

10. Perform statistical analysis of a given marker (eg PD1) MFI, for the different meta-clusters:
```R
StatAnalysisPD1 <- CytoSOM::BoxPlotMarkerMetaClust (CytoTree,Title = "MyTitle",tableTreatmentFCS,ControlTreatmen="PBS",BottomMargin=3,"PD1",Robust = TRUE,ClustHeat=TRUE)
```
A file `MyTitle_BoxPlotPD1MetaClust.pdf` is produced. Two files containing p-values are also produced: `MyTitle_LmPvalPD1Metacl.csv` and `MyTitle_PairwisePvalPD1Metacl.csv`
