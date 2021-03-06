% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FindMetaClustNames.R
\name{FindMetaClustNames}
\alias{FindMetaClustNames}
\title{Find each name of metaclusters that contains all sub-names}
\usage{
FindMetaClustNames(subNames, TreeMetaCl, start = F)
}
\arguments{
\item{subNames}{list of sub-part of a meta-cluster name}

\item{TreeMetaCl}{FlowSOM tree with meta-clusters, constructed within buildFSOMTree (from the data of FSOMData)}

\item{start}{If TRUE, match subnames with the start of meta-cluster names}
}
\value{
a vector of names
}
\description{
Find each name of metaclusters that contains all sub-names
}
