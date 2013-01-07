\name{rlme-package}
\alias{rlme-package}
\docType{package}
\title{
rlme
}
\description{
An R package for rank-based
estimation and prediction in random
effects nested models
}
\details{
\tabular{ll}{
Package: \tab rlme\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2013-01-07\cr
License: \tab GPL (>= 2)\cr
}

}
\author{
Yusuf Bilgic \email{bilgic@geneseo.edu} and Herb Susmann \email{hps1@geneseo.edu}

Maintainer: Yusuf Bilgic \email{bilgic@geneseo.edu}
}
\keyword{ models }
\keyword{ package }
\seealso{
\code{\link{rlme}}
}
\examples{

library(rlme)
data(schools)
rlme.fit = rlme(y ~ 1 + sex + age + (1 | school) + (1 | school:section), schools)
summary(rlme.fit)

}
