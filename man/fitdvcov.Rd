\name{fitdvcov}
\alias{fitdvcov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitdvcov
}
\description{
Obtains measurement for the fits based on estimates beta1, beta2 and covariance matrix from a rank based methods.
}
\usage{
fitdvcov(x1, beta1, beta2, vcw)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
%%     ~~Describe \code{x1} here~~
}
  \item{beta1}{
%%     ~~Describe \code{beta1} here~~
}
  \item{beta2}{
%%     ~~Describe \code{beta2} here~~
}
  \item{vcw}{
%%     ~~Describe \code{vcw} here~~
}
}
\seealso{
\code{\link{compare.fits}}
}
\examples{

# Compare GR and JR methods

data(schools)

model = y ~ 1 + sex + age + (1 | region) + (1 | region:school)

# Extract covariants into matrix
cov = as.matrix(data.frame(schools[,"sex"], schools[,"age"]))

# Fit the models using each method
jr.fit = rlme(model, schools, method="jr")
gr.fit = rlme(model, schools, method="gr")

# Extract beta estimates, ignoring the intercept
jr.beta = jr.fit$fixed.effects$Estimate[c(2, 3)]
gr.beta = gr.fit$fixed.effects$Estimate[c(2, 3)]

# Extract beta variance matrix
var.b = jr.fit$var.b

fitdvcov(cov, jr.beta, gr.beta, var.b)

}
