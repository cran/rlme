\name{JR_est}
\alias{JR_est}
\title{
JR Method
}
\description{
Fit a model using the JR method
}
\usage{
JR_est(x, y, I, sec, mat, school, section, rprpair = 'hl-disp')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Covariate matrix or data frame
}
  \item{y}{
    Response matrix or data frame
}
  \item{I}{
    Number of clusters.

%%     ~~Describe \code{I} here~~
}
  \item{sec}{
    A vector of subcluster numbers in clusters.

}
  \item{mat}{
  A matrix of numbers of observations in subclusters. Dimension is Ixmax(number ofsubclusters). Each row indicates one cluster. 
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
  \item{rprpair}{
  By default, it uses "hl-disp" in the random prediction procedure (RPP). Also, "med-mad" would be an alternative.
%%     ~~Describe \code{section} here~~
}

}
\value{
    
  \item{theta}{
    Fixed effect estimates.
}
  \item{ses}{
    Standard error for the fixed esimates.
}
  \item{sigma}{
    Covariate variance estimates using RPP (Groggel and Dubnicka's procedure).
}
  \item{ehat}{
    Raw error.
}
  \item{effect_sch}{
    Cluster random error.
}
  \item{effect_sec}{
    Subcluster random error.
}
  \item{effect_err}{
    Epsilon error.
}

}
\author{
Yusuf Bilgic
}
\seealso{
\code{\link{rlme}}
}
