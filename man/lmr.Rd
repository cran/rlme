\name{lmr}
\alias{lmr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rank Based Fixed Effect Regression
}
\description{
Computes rank based regression estimates for fixed effect models.
}
\usage{
lmr(f, data, se=FALSE, method='L-BFGS-B')
}
\arguments{
  \item{f}{
    A model formula
}
  \item{data}{
   Data to use for model fitting
}
  \item{se}{
    Boolean indicating whether or not to calculate standard errors for intercept and slope estimates
}
  \item{method}{
    Optimization method to use. Will accept any method usable by optim, e.g. one of c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
    "Brent"). "BFGS" or "L-BFGS-B" are reccomended. "L-BFGS-B" should be used for large datasets to conserve memory. 
}
}
\value{
  \item{fixed.effects}{
    Fixed effect estimates
}
  \item{ehat}{
    Residuals from model
}
}
\author{
Herb Susmann
}
\seealso{
\code{\link{rlme}}
\code{\link{optim}}
}
\examples{
# load schools data
data(schools)

# Fit fixed effects model with lmr
lmr.fit = lmr(y ~ age + sex, data=schools)

summary(lmr.fit)

# Fit with lmr and calculate standard errors
lmr.fit = lmr(y ~ age + sex, data=schools, se=TRUE)

summary(lmr.fit)
}

