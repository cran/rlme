\name{compare.fits}
\alias{compare.fits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compare Fits
}
\description{
Compares two model fits. It returns tdbeta value and cfits values of two fits. The function uses the fixed effects estimates from fit 1 and fit 2 along with the covariance of the rank-based fit. 
}
\usage{
compare.fits(x, fit1, fit2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Matrix of covariates
}
  \item{fit1}{
  A class of type rlme.
}
  \item{fit2}{
  A class of type rlme. 
}
}
\value{
Returns tdbeta and cfits values.
}
\seealso{
\code{\link{fitdvcov}}
}
\examples{

data(schools)
model = y ~ 1 + sex + age + (1 | region) + (1 | region:school)

# Extract covariants into matrix
cov = as.matrix(data.frame(schools[,"sex"], schools[,"age"]))

# Fit the models using each method
reml.fit = rlme(model, schools, method="reml")
gr.fit = rlme(model, schools, method="gr")

compare.fits(cov, reml.fit, gr.fit)
}
\keyword{ models }
