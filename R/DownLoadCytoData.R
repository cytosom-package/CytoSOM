#' Download fcs files and select the population of interest in FlowJo workspace
#' 
#' DownLoadCytoData is a tool to download data contained in fcs files and select the ones that have been gated as the population of interest in FlowJo. Folder where FCS files are localised and FlowJo workspace should be in current environment.
#' @param dirFCS sub-directory of working directory, where fcs files are located
#' @param gatingName name of the population of interest selected through a FlowJo gating
#' @param fcsPattern name of a pattern common to the fcs files (optional but useful to select only a selection of FCS files contained in a same folder)
#' @param compensate TRUE if data needs to be compensated
#' @return FlowSOM data object that includes the data contained in the FlowJo gating and the name of the population of interest.
#' @examples CytoData=DownLoadCytoData(dirFCS="Experiment1", gatingName="CD3+", fcsPattern="exp1", Compensate=FALSE)
#' @export
DownLoadCytoData <- function(dirFCS="",gatingName,fcsPattern = "Tube",compensate=FALSE ){
    flowJoWS=list.files(pattern=".wsp")
    if (length(flowJoWS) > 1)
    {
        stop("several FlowJoWorkSpace")
    }
    if (dirFCS == "") {absoluteDirFCS = getwd() }
    else {
        absoluteDirFCS=paste(getwd(),"/",dirFCS,"/",sep="")
        }
    files<-list.files(path = absoluteDirFCS , pattern=fcsPattern)
    if ((unlist(packageVersion("CytoML"))[1] == 1) & (unlist(packageVersion("CytoML"))[2] >= 12)) {
     ## print("Parse flowJo within version 1.12 of CytoML")
        data<-parse_flowjo_CytoML_v12(files,flowJoWS)}
    else {
            data<-parse_flowjo_CytoML(files,flowJoWS)
        }
    dataGated<-gating_subset_toolBox(data,gatingName)
    fSOM<-FlowSOM::ReadInput(dataGated$flowSet,compensate = compensate,transform = FALSE,scale = FALSE,scaled.center = TRUE,scaled.scale = TRUE,silent = FALSE)
    return(list(fSOMData=fSOM,flJoDataGated=dataGated,gatingName=gatingName))
}
