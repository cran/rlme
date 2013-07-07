\name{instruction}
\alias{instruction}
\docType{data}
\title{
Instruction
}
\description{
A data frame on school instruction results.
}
\usage{data(instruction)}
\format{
  A data frame with 1190 observations on the following 13 variables.
  \describe{
    \item{\code{X}}{a numeric vector}
    \item{\code{girl}}{a numeric vector}
    \item{\code{minority}}{a numeric vector}
    \item{\code{mathkind}}{a numeric vector}
    \item{\code{mathgain}}{a numeric vector}
    \item{\code{ses}}{a numeric vector}
    \item{\code{yearstea}}{a numeric vector}
    \item{\code{mathknow}}{a numeric vector}
    \item{\code{housepov}}{a numeric vector}
    \item{\code{mathprep}}{a numeric vector}
    \item{\code{classid}}{a numeric vector identifying the class within school}
    \item{\code{schoolid}}{a numeric vector identifying the school}
    \item{\code{childid}}{a numeric vector}
  }
}
\source{
West, B., Welch, K. B., & Galecki, A. T. (2006). Linear mixed models: a practical guide using statistical software. Chapman & Hall/CRC.
}
\examples{
# The following code takes a few minutes to run.
# In the interest of saving CRAN's example testing time,
# it has been commented out. If you want to use it,
# just uncomment and run.

# data(instruction)
# attach(instruction)

# data = data.frame(
#   y = mathgain,
#   mathkind = mathkind, 
#   girl = girl,
#   minority = minority,
#   ses = ses, 
#   school = factor(schoolid), 
#   section = factor(classid))


# fit.rlme = rlme(y ~ 1 + mathkind + girl + minority + ses + (1 | school) + (1 | school:section),
#  data = data,
#  method = "gr")
  
# summary(fit.rlme)

}
\keyword{datasets}
