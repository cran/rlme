\name{GEER_est}
\alias{GEER_est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GEER: General Estimating Equation Rank-Based Estimation Method
}
\description{
The package rlme calls this function for gee method, one of the methods proposed in Bilgic's study (2012). Also see Kloke et al. (2013). 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GEER_est(x, y, I, sec, mat, school, section, weight, rprpair, verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Design matrix, pxn, without intercept. 
%%     ~~Describe \code{x} here~~
}
  \item{y}{
  Response vector of nx1.
%%     ~~Describe \code{y} here~~
}
  \item{I}{
  Number of clusters.
%%     ~~Describe \code{I} here~~
}
  \item{sec}{
  A vector of subcluster numbers in clusters.
%%     ~~Describe \code{sec} here~~
}
  \item{mat}{
  A matrix of numbers of observations in subclusters.
  Dimension is Ixmax(number ofsubclusters). Each row indicates one cluster. 
%%     ~~Describe \code{mat} here~~
}
  \item{school}{
  A vector of clusters, nx1. 
%%     ~~Describe \code{school} here~~
}
  \item{section}{
  A vector of subclusters, nx1.
%%     ~~Describe \code{section} here~~
}

  \item{weight}{
When weight="hbr", it uses hbr weights in GEE weights. By default, ="wil", it uses Wilcoxon weights. See the theory in the references. 
%%     ~~Describe \code{section} here~~
}

  \item{rprpair}{
  By default, it uses "hl-disp" in the random prediction procedure (RPP). Also, "med-mad" would be an alternative.
%%     ~~Describe \code{section} here~~
}
  \item{verbose}{
    Boolean indicating whether to print out diagnostic messages.
}

}
\references{
Y. K. Bilgic. Rank-based estimation and prediction for mixed effects models in nested designs. 2012. URL http://scholarworks.wmich.edu/dissertations/40. Dissertation.

A. Abebe, J. W. McKean, J. D. Kloke and Y. K. Bilgic. Iterated reweighted rank-based estimates for gee
models. 2013. Submitted.

%% ~put references to the literature/web site here ~
}
\author{
Yusuf K. Bilgic, yekabe@hotmail.com
%%  ~~who you are~~
}

\seealso{
rlme, GR_est, JR_est, rprmeddisp
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# See the rlme function.
}
\keyword{ models }
