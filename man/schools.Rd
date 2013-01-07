\name{schools}
\alias{schools}
\docType{data}
\title{
PISA Literacy Data
}
\description{
The data in Program for International Assessment (PISA) on academic proficiency in schools around the world.
}
\usage{data(schools)}
\format{
  A data frame with 334 observations on the following 6 variables.
  \describe{
    \item{\code{y}}{a numeric vector indicating student literacy}
    \item{\code{socio}}{a numeric vector}
    \item{\code{sex}}{a numeric vector}
    \item{\code{age}}{a numeric vector}
    \item{\code{school}}{a numeric vector indicating the school}
    \item{\code{section}}{a numeric vector indicating the section within school}
  }
}
\references{
OECD (2010). PISA 2009 Results. http://www.oecd.org/
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
#
# The example takes a few seconds to run, so in order to save CRAN's testing time
# it has been commented out. To run, simply uncomment and execute.
#

# data(schools)
# rlme.fit = rlme(y ~ 1 + sex + age + (1 | school) + (1 | school:section), schools, method="gr")
# summary(rlme.fit)

}
\keyword{datasets}
