\name{rlme-package}
\alias{rlme-package}
\docType{package}
\title{
rlme
}
\description{
An R package for rank-based
robust estimation and prediction in random
effects nested models
}
\details{
\tabular{ll}{
Package: \tab rlme\cr
Type: \tab Package\cr
Version: \tab 0.2\cr
Date: \tab 2013-07-07\cr
License: \tab GPL (>= 2)\cr
}

}
\author{
Yusuf Bilgic \email{bilgic@geneseo.edu}, Herb Susmann \email{hps1@geneseo.edu} and Joseph McKean \email{joemckean@yahoo.com}

Maintainer: Yusuf Bilgic \email{bilgic@geneseo.edu} or \email{yusuf.k.bilgic@gmail.com}
}
\keyword{ models }
\keyword{ package }
\seealso{
\code{\link{rlme}}
}
\examples{

library(rlme)
data(schools)
formula = y ~ 1 + sex + age + (1 | region) + (1 | region:school)
rlme.fit = rlme(formula, schools)
summary(rlme.fit)
}
