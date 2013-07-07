\name{stanresidgr}
\alias{stanresidgr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Standard Residuals
}
\description{
Standardizes the residuals obtained from the GR fitting.
}
\usage{
stanresidgr(x, y, resid, delta = 0.8, param = 2, conf = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Design matrix.
%%     ~~Describe \code{x} here~~
}
  \item{y}{
  Response vector.
%%     ~~Describe \code{y} here~~
}
  \item{resid}{
  Residuals obtained from the rank-based fitting.
%%     ~~Describe \code{resid} here~~
}
  \item{delta}{
  See HM (2012).
%%     ~~Describe \code{delta} here~~
}
  \item{param}{
  See HM (2012).
%%     ~~Describe \code{param} here~~
}
  \item{conf}{
  See HM (2012).
%%     ~~Describe \code{conf} here~~
}
}
\references{
T. P. Hettmansperger and J. W. McKean. Robust Nonparametric Statistical Methods. Chapman Hall, 2012.

Y. K. Bilgic. Rank-based estimation and prediction for mixed effects models in nested designs. 2012. URL http://scholarworks.wmich.edu/dissertations/40. Dissertation.

%% ~put references to the literature/web site here ~
}
\author{
J. W.  McKean
}
