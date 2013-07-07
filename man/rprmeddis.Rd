\name{rprmeddis}
\alias{rprmeddis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rprmeddis
}
\description{
Robust rank-based prediction algorithm that gets predictions for random errors in three-level nested design. It needs one location and scale estimators. Hodges-Lehmann location estimate and dispersion functional estimate pair is called with rprpair="hl-disp" -by default- ; median and MAD pair is called with rprpair="med-mad" in rlme().

The rprmeddisp() function yields predictions of random effects and errors vectors along with scale estimates in each level. This function was designed for three-level nested design. See rprmeddisp2() in the package, this is for two-level nested design.

}
\usage{
rprmeddis(I, sec, mat, ehat, location, scale, rprpair = 'hl-disp')
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
  \item{ehat}{
  The residuals that inherits random effects and error effect to be predicted. 
%%     ~~Describe \code{ehat} here~~
}
  \item{location}{
    If location = scale = 1 then use Median and MAD in RPP
    If location = scale = 2 then use HL & Dispvar in RPP
    Note: this is deprecated. You should specify the location & scale parameters by using the rprpair parameter.
}
  \item{scale}{
  1 means mad, 2 means disp as scale estimators
}
  \item{rprpair}{
Character string indicating the location and scale parameters to use. Default to "hl-disp", but may also be "med-mad". See Bilgic (2012). 

}
}
\references{
Y. K. Bilgic. Rank-based estimation and prediction for mixed effects models in nested designs. 2012. URL http://scholarworks.wmich.edu/dissertations/40. Dissertation.

%% ~put references to the literature/web site here ~
}
\author{
Yusuf Bilgic \email{yekabe@hotmail.com}
}

\seealso{
\code{\link{rpr}}
\code{\link{dispvar}}
}
