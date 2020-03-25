
## User tool: Download data, given fcs files, FlowJo workspace should in in current environment, FCS directory is inside wd, given with no "/"
#' Download data, given fcs and FlowJo workspace
#' @param dirFCS sub-directory of working directory, where fcs files are
#' @param gatingName name of FlowJo gating to be applied
#' @param fcsPattern pattern of all fcs file names
#' @param compensate TRUE if data needs to be compensated?
#' @return FlowSOM data object
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
