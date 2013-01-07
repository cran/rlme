\name{LM_est}
\alias{LM_est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear Model Estimation using the nlme package.
}
\description{
This gets the REML or ML estimates and predictions of random effects from the nlme package.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LM_est(x, y, dat, method = "REML")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
   Design matrix, (p+1)xn, with intercept. 
%%     ~~Describe \code{x} here~~
}
  \item{y}{
  Response vector of nx1.
%%     ~~Describe \code{y} here~~
}
  \item{dat}{
  Data frame
}
  \item{method}{
  Character string indicating method to use, either "ML" or "REML" (defaults to REML).
}
}
\value{
  \item{theta}{
    Fixed effects esimates.
} 
  \item{ses}{
    Standard error for fixed effects.
} 
  \item{varb}{
  Variances.
}
  \item{sigma}{
  Error.
}
  \item{ehat}{
  Raw residuals
}
  \item{standr.lme}{
  Standardized residual
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

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
J. Pinheiro, D. Bates, S. DebRoy, D. Sarkar and R Development Core Team. nlme linear and non-
linear mixed effects models. The R Journal, 2011. URL http://CRAN.R-project.org/package=nlme.
R package version 3.1-98.
%% ~put references to the literature/web site here ~
}
\author{
Yusuf Bilgic
}
\seealso{
\code{\link{rlme}}
}
\keyword{ models }
