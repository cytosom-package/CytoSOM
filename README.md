CytoSOM: a package for easy use of FlowSOM
==========================================
Introduction
==============
CytoSOM helps to apply [FlowSOM](https://github.com/SofieVG/FlowSOM) to flow cytometry data, and perform statistical analysis. You can use the [userguide](https://github.com/gautierstoll/CytoSOM/blob/master/CytoSOM%20User%20Guide.pdf) or read the "HowTo" below.

## Contact
cytosom.package@gmail.com

## Installation

```R
devtools::install_github(repo ="gautierstoll/CytoSOM")
```

## HowTo

`CytoSOM` is applied when data has been gated within FlowJo software. A gating can be construct manually, see below.

`CytoSOM` needs to be installed. In a folder (eg `MyData`), put your FlowJo sesssion file (eg `file.wsp`) and a sub-folder (eg `FCSdata`) that contains the `.fcs` files (eg `Tube`...`.fcs`)

1. Launch `RStudio`

2. Set you working direcory as the folder where your FlowJo session and you data are, using the full folder name:
```R
setwd("full_name_before_myData/MyData")
```

3. Create an R object that collect your data (eg `CytoData`), indicating the name of the gating used in FlowJo (eg `CD45`):
```R
CytoData <- CytoSOM::DownLoadCytoData(dirFCS="FCSdata","CD45",fcsPattern = "Tube",compensate=FALSE)
```
Note that in this case the data have been already compensated.

4. Build a clustering tree (eg `CytoTree`), indicating the list of makers used for clustering (eg `c("CD4","CD8","CD11b","FOXP3")`), the size of the cluster grid (eg `7`), the number of meta-clusters (eg `10`), and the seed of the random generator (eg `0`):
```R
CytoTree <- CytoSOM::buildFSOMTree(CytoData,c("CD3","CD4","CD8","CD11b","FOXP3", "CD19"),7,25,0)
```

5. Plot the tree
```R
CytoSOM::PlotStarsMSTRm(CytoTree$fSOMTree,CytoTree$metaCl,"Title Name",0)
```
If the smallest clusters make the image diffcult to interpret, you can plot the tree by removing the smallest clusters (eg remove the 3 smallest clusters)
```R
CytoSOM::PlotStarsMSTRm(CytoTree$fSOMTree,CytoTree$metaCl,"Title Name",3)
```
If the tree look staitfactory, you can do the statsistical analysis. Otherwise, try to rebuild the tree with other paramters (cluster grid size and/or meta-cluster number)

6. Rename the meta-cluster names, using a set of markers (eg `c("CD4","CD8","CD11b","FOXP3", "CD19")`)
```R
CytoTree <- CytoSOM::TreeMetaRenaming(CytoTree,c("CD4","CD8","CD11b","FOXP3", "CD19"),"shortRobustName")
```

7. Download an annotation table in `.csv` format (eg `AnnotationTable.csv`), that contains a column 'files', a column 'Treatment' and an optinal column 'NormalizationFactor', indicating the separator (eg `;`):
```R
tableTreatmentFCS <- read.csv("AnnotationTable.csv",sep=";")
```

8. Plot trees in `.pdf` files, for different population sizes and different markers (eg `c("MHCII","PD1","PDL1","PDL2")`):
```R
CytoSOM::plotTreeSet(CytoTree ,c("MHCII","PD1","PDL1","PDL2"),"TitleTrees",rmClNb=0,tableTreatmentFCS,globalScale=T)
```
This produces two `.pdf` files in the working directory: `TitleTrees_ClusterTree.pdf` and `TitleTrees_MarkerTree.pdf`.

9. Perform statistical analysis of meta-cluster sizes:
```R
StatAnalysisSizes <- CytoSOM::BoxPlotMetaClust(CytoTree,Title="MyTitle",tableTreatmentFCS,ControlTreatmen="PBS",
BottomMargin=3,yLab="CD45",Norm=FALSE,Robust = TRUE,ClustHeat=TRUE)
```
If `Norm` is set to `TRUE`, the column 'NormalizationFactor' is used to normalize the meta-cluster sizes. Otherwise, the analysis is perfomed on relative size (precentage). A file `MyTitle_BoxPlotPercentMetaClust.pdf` is produced. Two files containing p-values are also produced: `MyTitle_LmPvalPercentMetacl.csv` and `MyTitle_PairwisePvalPercentMetacl.csv`

10. Perform statistical analysis of a given marker (eg PD1) MFI, for the different meta-clusters:
```R
StatAnalysisPD1 <- CytoSOM::BoxPlotMarkerMetaClust(CytoTree,Title="MyTitle",tableTreatmentFCS,ControlTreatmen="PBS",
BottomMargin=3,"PD1",Robust = TRUE,ClustHeat=TRUE)
```
A file `MyTitle_BoxPlotPD1MetaClust.pdf` is produced. Two files containing p-values are also produced: `MyTitle_LmPvalPD1Metacl.csv` and `MyTitle_PairwisePvalPD1Metacl.csv`

### Gating witout FlowJo

`CytoSOM` needs to be installed. In a folder (eg `MyData`), put a sub-folder (eg `FCSdata`) that contains the `.fcs` files (eg `Tube`...`.fcs`)

1. Launch `RStudio`

2. Set you working directory as the folder where you data are, using the full folder name:
```R
setwd("full_name_before_myData/MyData")
```

3. Download the data:

```R
RawData <- FlowSOM::ReadInput(input = "FCSdata",pattern = "Tube",compensate = F)
```
Note that in this case the data have been already compensated.

4. Create two polygon gates, eg the first one within the 2D plot "FSC-A" x "SSC-A" using `.fcs` files 1 and 3, the second one within the 2D plot "FSC-A" x "Livedead" (from 0 to 10000) using `.fcs` files 2 and 4:
```R
Poly1 <- CytoSOM::InteractivePolyGate(RawData,marker1 = "FSC-A",marker2 = "SSC-A",fcsFiles = c(1,3))
Poly2 <- CytoSOM::InteractivePolyGate(RawData,marker1 = "FSC-A",marker2 = "Livedead",fcsFiles = c(2,4),ylim=c(0,10000))
```

5. Create a dataset with the instersection of the two gates above (named "CD45"):
```R
CytoData <- CytoSOM::PolygonGatingRawData(RawData,Polygons = list(Poly1,Poly2),gatingName = "CD45”)
```

6. A new polygon can be applied to this gated dataset:

```R
Poly3 <- CytoSOM::InteractivePolyGate(CytoData$fSOMData,marker1 = "SSC-H",marker2 = "Livedead",fcsFiles = c(8,10))
CytoData <- CytoSOM::PolygonGatingGatedData(CytoData,Polygons = list(Poly3),gatingName = "CD45”)
```

Then the analysis can be continued at point 4 above.
