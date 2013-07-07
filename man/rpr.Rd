\name{rpr}
\alias{rpr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Frontend to rprmeddis
}
\description{
Used as a simple wrapper for rprmeddis
}
\usage{
rpr(ehat, school, section, rprpair = "hl-disp")
}
\arguments{
  \item{ehat}{
    The raw residual. 

%%     ~~Describe \code{ehat} here~~
}
  \item{school}{
    The number of clusters.

%%     ~~Describe \code{school} here~~
}
  \item{section}{
    The number of subclusters.

%%     ~~Describe \code{section} here~~
}
  \item{rprpair}{
Character string indicating the location and scale parameters to use. Default to "hl-disp", but may also be "med-mad". See Bilgic (2012). 
}
}
\references{
Y. K. Bilgic. Rank-based estimation and prediction for mixed effects models in nested designs. 2012. URL http://scholarworks.wmich.edu/dissertations/40. Dissertation.
}
\author{
J. W. McKean and Y. K. Bilgic
}
\seealso{
\code{\link{rprmeddis}}
\code{\link{dispvar}}
}
