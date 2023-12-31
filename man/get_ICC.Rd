% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/commonCorrelationModel.R
\name{get_ICC}
\alias{get_ICC}
\title{Intraclass correlation}
\usage{
get_ICC(n_neither, n_one, n_both, method = "auto", print.out = TRUE)
}
\arguments{
\item{n_neither}{a numeric, number of samples in which neither of the twins showed the relevant phenotype.}

\item{n_one}{a numeric, number of samples in which only one of the twins showed the relevant phenotype.}

\item{n_both}{a numeric, number of samples in which both twins showed the relevant phenotype.}

\item{method}{character, methods for parameter estimation (e.g. "auto","normal","gof"). Default is "auto".
Recommended method is 'normal' when total number of sample N >100 and 
the prevalence parameter pi (0.3<pi<0.7) and the level of twin similarity
with respect to the dichotomous trait in question Rho is no greater than 0.60.
Recommended method is 'gof' when total number of sample N is small or moderate in size (<=100), 
pi is fairly extreme (<0.3 or >0.7), or p is anticipated to be of substantial magnitude (>=0.6)}

\item{print.out}{a boolean, Wwether or not to print results.}
}
\value{
Return a dataframe. 
pi is the estimated prevalence of trait.
Vh is the standard error square of ICC.
ICC means the the level of twin similarity with respect to the dichotomous trait.
Z0 and p_value is the statistic using the normal theory method.
CI_lower and CI_upper is the 95% confidence interval of ICC.
}
\description{
\code{get_ICC} is a function measures twin similarity with the intraclass correlation (ICC) 
of the phenotype separately within monochorionic and dichorionic pairs.
}
\examples{
data <- load_example_data();
chisq_res <- twin.chisqTest(data = data,famid = 'FAMID',outcome = 'BPD',zyg = 'zygo')
stat_mat <- chisq_res$stat_mat
get_ICC(n_neither = stat_mat['MZ','Neither'],
        n_one     = stat_mat['MZ','One'],
        n_both    = stat_mat['MZ','Both'])

}
