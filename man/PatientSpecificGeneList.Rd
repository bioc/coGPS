\name{PatientSpecificGeneList}
\alias{PatientSpecificGeneList}
\title{Patient Specific outlier gene list}
\description{Generate an outlier gene list for each patient restricted to the top PCOPA scored genes}
\usage{
PatientSpecificGeneList(exprslist, alpha, side, type,TopGeneNum)
}
\arguments{
  \item{exprslist}{Each element of  \eqn{exprslist} is a list with the first element being  \eqn{exprs} and the second element being  \eqn{classlab}. Each row of  \eqn{exprs} represents one gene and each column represents one sample.  \eqn{classlab} is a zero-one vector indicating the status of samples. We use 0 for the baseline group, usually the normal group, and 1 for the comparison group, usually the tumor group.}
  \item{alpha}{Significance level for P-value.}
  \item{side}{ A vector specifying the definition of P-value in each of the study, which could be either  \eqn{up},  \eqn{down}, or  \eqn{twosided}.}
  \item{type}{ A vector specifying whether the outlier pattern is  \eqn{subtype} or  \eqn{uniform}.}
  \item{TopGeneNum}{ a number specifying the top number of outlier genes scored by PCOPA to be included in the generation of individual outlier gene list for each patient.}
}

\value{
\item{outliergene_bypatient}{ a list whose length equals the number of tumor samples (patients). each element of the list is a list of length equaling to the length of \eqn{exprslist}, in other words the number of studies(or data type), showing the outlier gene for each patient in each study (or data type)}
}
\references{
Wei, Y., Hennessey, P., Gaykalova, D., Califano, J.A., Ochs, M.F., (2011) Cancer Outlier Gene Profile Sets Elucidate Pathways in Head and Neck Squamous Cell Carcinoma.
}
\author{Yingying Wei}

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

#generate an outlier gene list for each patient restricted to the top PCOPA scored genes
IndividualList7<-PatientSpecificGeneList(trylist,0.05,side=c("up","down","up"),type="subtype",TopGeneNum=100)
}
