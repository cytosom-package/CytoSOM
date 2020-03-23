MetaClusterNaming <- function(TreeMetaCl,Markers)
{
    MarkerIn = sapply(Markers,function(Marker){length(which(TreeMetaCl$fSOMTree$prettyColnames == Marker))})
    Markers = Markers[which(MarkerIn > 0)]
    print(paste("Use Marker",Markers))


    metaClustListName=lapply(unique(TreeMetaCl$metaCl),function(metaClust){
            clusterList=which(TreeMetaCl$metaCl == metaClust)
            metaClustIndices=unlist(sapply(clusterList,function(cluster){which(TreeMetaCl$fSOMTree$map$mapping[,1] == cluster)}))
            nameList = sapply(Markers,function(Marker){
                MarkerIndex=which(TreeMetaCl$fSOMTree$prettyColnames == Marker)
                metaClustMedian=median(TreeMetaCl$fSOMTree$data[metaClustIndices,MarkerIndex],na.rm=T)
                vectQuantiles=quantile(TreeMetaCl$fSOMTree$data[,MarkerIndex],na.rm=T)
                nonRobustName = c(paste(Marker,"-",sep=""),paste(Marker,"+",sep=""))[as.numeric(metaClustMedian > vectQuantiles[3])+1]
                robustName = c(paste(Marker,"-",sep=""),paste(Marker,"med",sep=""),paste(Marker,"+",sep=""))[as.numeric(metaClustMedian > vectQuantiles[2])+as.numeric(metaClustMedian > vectQuantiles[4])+1]
                robustNameNoMed= c(paste(Marker,"-",sep=""),"",paste(Marker,"+",sep=""))[as.numeric(metaClustMedian > vectQuantiles[2])+as.numeric(metaClustMedian > vectQuantiles[4])+1]
                return(c(nonRobustName,robustName,robustNameNoMed))
            })
            c(metaClust,apply(nameList,1,function(l){paste(l,collapse="")}))
    })
    metaClustDF=as.data.frame(do.call(rbind,metaClustListName),stringsAsFactors=F)
    names(metaClustDF)=c("number","nonRobustName","robustName","shortRobustName")
    return(metaClustDF)
}
