#' Plot a set of trees, in a pdf file, for metaclusters and markers
#'
#' plotTreeSet is a function that allows the automatic export of
#' \itemize{
#'  \item{one PDF file containing the plot of the general tree and its representation according to each experimental condition}
#'  \item{one PDF file containing the representation on the general and treatment specific of the median fluorescent intensity of one choosen marker}
#' }
#' the treatmentTable should be a dataframe with two column: "Treatment", "files"
#' @param TreeMetacl FlowSOM tree with meta-clusters, constructed by buildFSOMTree
#' @param markers vector of markers for which PDF files of MFI will be generated
#' @param Title title
#' @param rmClNb number of smallest clusters to remove
#' @param treatmentTable data frame containing a column 'files', a column 'Treatment' (case sensitive) and a column "NormaizationFactor"
#' @param globalScale TRUE if color scale for markers is the same for each treatments
#' @examples plotTreeSet(TreeMetaCl=CytoTree,markers=c("PD-1","CTLA-4"),Title="Experiment1",rmClNb=0,treatmentTable="Experiment1Table.csv",globalScale=TRUE)
#' @export
plotTreeSet <- function(TreeMetacl,markers,Title,rmClNb,treatmentTable,globalScale=T){
    if (rmClNb>0) {
        indexKeep =  which(TreeMetacl$fSOMTree$MST$size > sort(TreeMetacl$fSOMTree$MST$size)[rmClNb])} else {indexKeep = 1:length(TreeMetacl$fSOMTree$MST$size)}
    pdf(file=NoSpCharForFile(paste(Title,"_ClusterTree.pdf",sep="")))
    ## plot tree of pooled data
    PlotStarsMSTRm(TreeMetacl$fSOMTree,TreeMetacl$metaCl,Title,rmClNb)
    PlotStarsMSTRm(TreeMetacl$fSOMTree,TreeMetacl$metaCl,paste(Title," (id. sizes)",sep=""),rmClNb,equalSize=T)
    PlotStarsMSTRm(TreeMetacl$fSOMTree,TreeMetacl$metaCl,paste(Title," (large legend)",sep=""),rmClNb,smallTree=T)
    if (is.null(TreeMetacl$metaClNumber)) {
        PlotLabelsRm(TreeMetacl$fSOMTree,TreeMetacl$metaCl,paste(Title,"\nMetaclusterTree",sep=""),rmClNb,equalSize=T)}
    else {PlotLabelsRm(TreeMetacl$fSOMTree,TreeMetacl$metaClNumber,paste(Title,"\nMetaclusterTree",sep=""),rmClNb,equalSize=T)}

    Treatments=unique(treatmentTable$Treatment[which(sapply(treatmentTable$files,function(files){length(grep(files,names(TreeMetacl$fSOMTree$metaData),fixed=T))>0}))])
    print("Treatments:")
    print(Treatments)

    ## plot tree of subdata for each treatments
    for (treatName in Treatments) {
        fcsFiles=treatmentTable$files[which(treatmentTable$Treatment == treatName)]
        treatIndex = unlist(sapply(fcsFiles,function(file){grep(file,names(TreeMetacl$fSOMTree$metaData),fixed=T)}))
        print(paste("Plot tree for treatment",treatName))
        PlotStarsMSTCondRm(TreeMetacl$fSOMTree,TreeMetacl$metaCl,treatIndex,paste(Title," Treat: ",treatName,sep=""),rmClNb)
    }
    dev.off()
    pdf(file=NoSpCharForFile(paste(Title,"_MarkerTree.pdf",sep="")))
    markersInData = intersect(markers,TreeMetacl$fSOMTree$prettyColnames)
    print("Markers:")
    print(markersInData)
    ##plot tree for each marker
    for (marker in markersInData)
    {
        uglyName = names(which(TreeMetacl$fSOMTree$prettyColnames == marker))
        if (globalScale) ## construct the global scale for color scale of markers
        {
            listTreatmentIndices = lapply(unique(treatmentTable$Treatment),function(treatName){
                fcsFiles=treatmentTable$files[which(treatmentTable$Treatment == treatName)]
                CondIndex = unlist(sapply(fcsFiles,function(file){grep(file,names(TreeMetacl$fSOMTree$metaData),fixed=T)}))
                return(unlist(sapply(TreeMetacl$fSOMTree$metaData[CondIndex],function(x){x[1]:x[2]})))
            })
            maxMin = matrix(unlist(lapply((1:length(TreeMetacl$fSOMTree$map$medianValues[,1]))[indexKeep],function(cluster){
                medianList=sapply(listTreatmentIndices,function(indices){
                    median(TreeMetacl$fSOMTree$data[intersect(indices,which(TreeMetacl$fSOMTree$map$mapping[,1] == cluster)),uglyName])});
                c(max(medianList,na.rm=T),min(medianList,na.rm=T))
            })),nrow=2)
            maxGlobal = max(maxMin[1,])
            minGlobal = min(maxMin[2,])
            print(paste("Scale for marker",marker,":",minGlobal,"->",maxGlobal))
            globMinMax=c(minGlobal,maxGlobal)
        }
        else {globMinMax=c()}
        ##print(paste("Marker: ",marker,sep=""))
        print(paste("Plot tree for marker",marker," -- ",uglyName))
        PlotMarkerMSTRm(TreeMetacl$fSOMTree,uglyName,paste(Title," Marker: ",marker,sep=""),rmClNb,globMinMax)
    }
    for (marker in markersInData){ ## plot tree for each marher and each treatment
        uglyName = names(which(TreeMetacl$fSOMTree$prettyColnames == marker))
        if (globalScale) ## construct the global scale for color scale of markers
        {
            listTreatmentIndices = lapply(unique(treatmentTable$Treatment),function(treatName){
                fcsFiles=treatmentTable$files[which(treatmentTable$Treatment == treatName)]
                CondIndex = unlist(sapply(fcsFiles,function(file){grep(file,names(TreeMetacl$fSOMTree$metaData),fixed=T)}))
                return(unlist(sapply(TreeMetacl$fSOMTree$metaData[CondIndex],function(x){x[1]:x[2]})))
            })
            maxMin = matrix(unlist(lapply((1:length(TreeMetacl$fSOMTree$map$medianValues[,1]))[indexKeep],function(cluster){medianList=sapply(listTreatmentIndices,function(indices){median(TreeMetacl$fSOMTree$data[intersect(indices,which(TreeMetacl$fSOMTree$map$mapping[,1] == cluster)),uglyName])});c(max(medianList,na.rm=T),min(medianList,na.rm=T))})),nrow=2)
            maxGlobal = max(maxMin[1,])
            minGlobal = min(maxMin[2,])
            globMinMax=c(minGlobal,maxGlobal)
        }
        else {globMinMax=c()}
        ##print(paste("Marker: ",marker,sep=""))
        for (treatName in Treatments){
            ##print(paste("Treatment: ",treatName,sep=""))
            fcsFiles=treatmentTable$files[which(treatmentTable$Treatment == treatName)]
            treatIndex = unlist(sapply(fcsFiles,function(file){grep(file,names(TreeMetacl$fSOMTree$metaData),fixed=T)}))
            print(paste("Plot tree for marker",marker," -- ",uglyName,"and treatment",treatName))
            PlotMarkerMSTCondRm(TreeMetacl$fSOMTree,uglyName,treatIndex,paste(Title," Marker: ",marker," Treat: ",treatName,sep=""),rmClNb,globMinMax)
        }
    }
    dev.off()
}
