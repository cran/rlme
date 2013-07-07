\name{wilstep}
\alias{wilstep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wilcoxon One Step Rank-based Estimate in GR Method
}
\description{
Gets weighted rank based fittings using the ww package. The Rfit package would be used as an alternative. Also, see the wilstep2 function in R codes for two-level nested designs.
}
\usage{
wilstep(I, sec, mat, init = F, y, x, sigmaa2 = 1,
        sigmaw2 = 1, sigmae2 = 1, thetaold = c(0),
        eps = 1e-04, iflag2 = 0, rprpair='hl-disp')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{I}{
    Number of clusters.

%%     ~~Describe \code{I} here~~
}
  \item{sec}{
    A vector of subcluster numbers in clusters.

%%     ~~Describe \code{sec} here~~
}
  \item{mat}{
   A matrix of numbers of observations in subclusters. Dimension is Ixmax(number ofsubclusters). Each row indicates one cluster. 
%%     ~~Describe \code{mat} here~~
}
  \item{init}{
%%     ~~Describe \code{init} here~~
}
  \item{y}{
    Response vector of nx1.

%%     ~~Describe \code{y} here~~
}
  \item{x}{
    Design matrix, pxn, without intercept. 

%%     ~~Describe \code{x} here~~
}
  \item{sigmaa2}{
  Initial sigma for cluster in three-level design.
%%     ~~Describe \code{sigmaa2} here~~
}
  \item{sigmaw2}{
    Initial sigma for subcluster in three-level design.

%%     ~~Describe \code{sigmaw2} here~~
}
  \item{sigmae2}{
    Initial sigma for error in three-level design.

%%     ~~Describe \code{sigmae2} here~~
}
  \item{thetaold}{
  Initial input.
%%     ~~Describe \code{thetaold} here~~
}
  \item{eps}{
  Epsilon value
%%     ~~Describe \code{eps} here~~
}
  \item{iflag2}{
  y or n
%%     ~~Describe \code{iflag2} here~~
}
  \item{rprpair}{
   Either 'hl-disp' or 'med-mad'
}
}
\details{
Initial inputs are from the independent model.
}
\references{
Y. K. Bilgic and J. W. McKean. Iteratively reweighted generalized rank-based method in mixed models. 2013. Under preperation.

J. T. Terpstra and J. W. McKean. Rank-based analysis of linear models using R. Journal of Statistical Software, 14(7) 1 - 26, 7 2005. ISSN 1548-7660. URL http://www.jstatsoft.org/v14/i07.

%% ~put references to the literature/web site here ~
}
\author{
J. W. McKean and Y. K. Bilgic
}
