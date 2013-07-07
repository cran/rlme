\name{beta_var}
\alias{beta_var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimate fixed-effect variance for Joint Rank Method (JR) in three-level nested design.
}
\description{
  Fixed effect variance estimation for Joint Rank Method (JR). It assumes Compound Symmetric (CS) structure of error terms. For k-level design, there are k-1 intra/inter-class parameters to place in a correlation matrix of errors. 
}
\usage{
beta_var(x, school, tauhat, v1, v2, v3, section, mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Data frame of covariates.
}
  \item{school}{
  A vector of cluster.
%%     ~~Describe \code{school} here~~
}
  \item{tauhat}{
  This is obtained from Rank-based fitting.
%%     ~~Describe \code{tauhat} here~~
}
  \item{v1}{
  This is 1, main diagonal element for correlation matrix of observations. Correlation of an observation with itself is 1.
%%     ~~Describe \code{v1} here~~
}
  \item{v2}{
  Intra-cluster correlation coefficient.
%%     ~~Describe \code{v2} here~~
}
  \item{v3}{
  Intra-subcluster correlation coefficient.
%%     ~~Describe \code{v3} here~~
}
  \item{section}{
    A vector of subclusters, nx1.
%%     ~~Describe \code{v3} here~~
}
  \item{mat}{
  A matrix of numbers of observations in subclusters. Dimension is Ixmax(number ofsubclusters). Each row indicates one cluster. 
%%     ~~Describe \code{mat} here~~
}
}
\details{
Correlation coefficients are obtained using Moment Estimates. See Klole et. al (2009), Bilgic (2012) and HM (2012)
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{var}{
    The variance of fixed estimated.
}
}
\references{
Y. K. Bilgic. Rank-based estimation and prediction for mixed effects models in nested designs.
2012. URL http://scholarworks.wmich.edu/dissertations/40. Dissertation.

J. Kloke, J. W. McKean and M. Rashid. Rank-based estimation and associated inferences for linear models with cluster correlated errors. Journal of the American Statistical Association, 104(485):384-390, 2009.

T. P. Hettmansperger and J. W. McKean. Robust Nonparametric Statistical Methods. Chapman Hall, 2012.
%% ~put references to the literature/web site here ~
}
\author{
Yusuf Bilgic
}
