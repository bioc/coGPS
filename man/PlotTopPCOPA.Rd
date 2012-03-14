\name{PlotTopPCOPA}
\alias{PlotTopPCOPA}
\title{Plot expression patterns of top ranked genes.}
\description{
It first sorts the expression value  \eqn{exprslist[[i]]\$exprs[j,]} among the baseline samples(e.g. normal ones) and comparison group (e.g. tumor ones)seperately for selected gene  \eqn{j}, and then plot the sorted expression values. 
The first argument  \eqn{exprslist} should be the same one as for \eqn{PCOPA}; the second argument \eqn{PCOPAresult} should be an output of PCOPA; the third argument  \eqn{topcut} determines how far we would go down the top ranked list; and the last argument  \eqn{typelist} is a vector specifying the titles for each graph corresponds to a specific study.
}
\usage{
PlotTopPCOPA(exprslist, PCOPAresult, topcut, typelist)
}

\arguments{
  \item{exprslist}{Each element of  \eqn{exprslist} is a list with the first element being  \eqn{exprs} and the second element being  \eqn{classlab}. Each row of  \eqn{exprs} represents one gene and each column represents one sample.  \eqn{classlab} is a zero-one vector indicating the status of samples. 
We use 0 for the baseline group, usually the normal group, and 1 for the comparison group, usually the tumor group.}
  \item{PCOPAresult}{Output of PCOPA.}
  \item{topcut}{ Cutoff of top ranked gene list.}
  \item{typelist}{ A vector specifying the titles for each graph corresponds to a specific study.}
}

\author{
Michael Ochs, Yingying Wei
}

\examples{
#read in data
data(Exon_exprs_matched)
data(Methy_exprs_matched)
data(CNV_exprs_matched)
data(Exon_classlab_matched)
data(Methy_classlab_matched)
data(CNV_classlab_matched)
head(Exon_exprs_matched)

#exprslist[[i]]$exprs should be in matrix format
Exon_exprs<-as.matrix(Exon_exprs_matched)
Methy_exprs<-as.matrix(Methy_exprs_matched)
CNV_exprs<-as.matrix(CNV_exprs_matched)

#exprslist[[i]]$classlab should be in vector format
Exon_classlab<-unlist(Exon_classlab_matched)
Methy_classlab<-unlist(Methy_classlab_matched)
CNV_classlab<-unlist(CNV_classlab_matched)

#make an exprslist consisting 3 studies
trylist<-list()
trylist[[1]]<-list(exprs=Exon_exprs,classlab=Exon_classlab)
trylist[[2]]<-list(exprs=Methy_exprs,classlab=Methy_classlab)
trylist[[3]]<-list(exprs=CNV_exprs,classlab=CNV_classlab)

#calculate P-value based statistics for outlier gene detection and output the outlier gene list for each patient
a7<-PCOPA(trylist,0.05,side=c("up","down","up"),type="subtype")

#plot expression patterns of top ranked genes. 
PlotTopPCOPA(trylist,a7,topcut=1,typelist=c("Exon","Methy","CNV"))
}