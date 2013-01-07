\name{plot.rlme}
\alias{plot.rlme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot rlme Fit
}
\description{
Generates Normal Q-Q plot of residuals from rlme fit
}
\usage{
\method{plot}{rlme}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  A list of class rlme. Store as fit.rlme.
}
  \item{\dots}{
  
}
}
\examples{
data(schools)
rlme.fit = rlme(y ~ 1 + sex + age + (1 | school) + (1 | school:section), schools, method="gr")
plot(rlme.fit)
}
