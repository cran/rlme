\name{hbrwts_gr}
\alias{hbrwts_gr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
HBR Weight
%%  ~~function to do ... ~~
}
\description{
Calculates hbr weights for the GEER method. This turns a vector of weights for a vector of errors. Used to make factor space more robust, up to 50\% breakdown. See HM (2012) and Terpstra and McKean (2005) for details. The ww package produces this weights as well.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hbrwts_gr(xmat, y, robdis2 = mycov.rob(as.matrix(xmat), method = "mcd")$robdis2,
          percent = 0.95, intest = myltsreg(xmat, y)$coef)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xmat}{
  Design matrix, pxn, without intercept. 
%%     ~~Describe \code{xmat} here~~
}
  \item{y}{
 Response vector in nx1.
%%     ~~Describe \code{y} here~~
}
  \item{robdis2}{
 This is obtained from mycov.rob(as.matrix(xmat, method = "mcd")$robdis2
%%     ~~Describe \code{robdis2} here~~
}
  \item{percent}{
 This is 0.95.
%%     ~~Describe \code{percent} here~~
}
  \item{intest}{
This is obtained from myltsreg(xmat, y)$coef
%%     ~~Describe \code{intest} here~~
}
}
\details{
The ww pack explains how it is obtained.
%%  ~~ If necessary, more details than the description above ~~
}
\references{
T. P. Hettmansperger and J. W. McKean. Robust Nonparametric Statistical Methods. Chapman Hall, 2012.

J. T. Terpstra and J. W. McKean. Rank-based analysis of linear models using R. Journal of Statistical
Software, 14(7):1 - 26, 7 2005. ISSN 1548-7660. URL http://www.jstatsoft.org/v14/i07.
%% ~put references to the literature/web site here ~
}
\author{
J. W. McKean
}
\seealso{
GEER_est
}
