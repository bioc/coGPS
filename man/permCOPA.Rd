\name{permCOPA}
\alias{permCOPA}
\title{ Calculate PCOPA value for permuations}
\description{ Run permutations by randomly shuffling the sample class labels and calculate a vector of PCOPA values for each permutation.
}
\usage{
permCOPA(exprslist, alpha=0.05, side, type, perms=100)
}
\arguments{
  \item{exprslist}{Each element of  \eqn{exprslist} is a list with the first element being  \eqn{exprs} and the second element being  \eqn{classlab}. Each row of  \eqn{exprs} represents one gene and each column represents one sample.  \eqn{classlab} is a zero-one vector indicating the status of samples. We use 0 for the baseline group, usually the normal group, and 1 for the comparison group, usually the tumor group.}
  \item{alpha}{Significance level for P-value.}
  \item{side}{ A vector specifying the definition of P-value in each of the study, which could be either  \eqn{up},  \eqn{down}, or  \eqn{twosided}.}
  \item{type}{ A vector specifying whether the outlier pattern is  \eqn{subtype} or  \eqn{uniform}.}
  \item{perms}{ Number of permutations to run.}
}

\value{
\item{permResult}{ A matrix where each row correspond to a gene and each column correspond to one permutation.}
}
\references{
Wei, Y., Hennessey, P., Gaykalova, D., Califano, J.A., Ochs, M.F., (2011) Cancer Outlier Gene Profile Sets Elucidate Pathways in Head and Neck Squamous Cell Carcinoma.
}
\author{
Michael Ochs
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

#run 2 permutations 
perma7<-permCOPA(trylist,0.05,side=c("up","down","up"),type="subtype",perms=2)
}
