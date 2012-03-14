\name{SampleData}
\alias{Exon_exprs_matched}
\alias{Exon_classlab_matched}
\alias{Methy_exprs_matched}
\alias{Methy_classlab_matched}
\alias{CNV_exprs_matched}
\alias{CNV_classlab_matched}
\alias{Hs.gmtl.c1}
\title{Sample Data for coGPS}
\description{Here we present an example of coGPS analysis.}

\arguments{
	\item{Exon_exprs_matched }{ Expression data for 44 tumors and 25 normals. Each row indicates a gene with row name showing gene name and each column indicates a sample with column name showing sample name.}
	\item{Exon_class_matched}{ A length 69 vector showing status of corresponding exon samples, 0 for normals and 1 for tumors.}
	\item{Methy_exprs_matched}{ Methylation data for 44 tumors and 25 normals.}
	\item{Methy_class_matched}{ A length 69 vector showing status of corresponding methylation samples, 0 for normals and 1 for tumors.}
	\item{CNV_exprs_matched}{ Copy number data for 44 tumors and 25 normals.}
	\item{CNV_class_matched}{ A length 69 vector showing status of corresponding copy number samples, 0 for normals and 1 for tumors.}
	\item{Hs.gmtl.c1}{ Broad Institute C1 Positional Gene Sets.}
}

\details{
In this application, the columns of each data type are matched. In other words, the first columns of Exon_exprs_matched, Methy_exprs_matched and CNV_exprs_matched correspond to the same patient. And hence the Exon_class_matched, Methy_class_matched and CNV_class_matched are identical. However, suppose in applications that we are not concerned with the outlier gene list for each patient, we can leave with the samples (columns) unmatched. 
}

