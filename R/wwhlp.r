#The ww package help files. See the authors.
#Not used in rlme
# Help file of R Functions for Weighted Wilcoxon Estimation and Inference
#---------------------------------------------------------------------------------------------
#Function:
#Description:
#Usage:
#Arguments:
#Value:
#Details:
#References:
#See Also:
#Examples:
#---------------------------------------------------------------------------------------------
#  
# Version:      1.6
# Last Update:  09/08/08
# 
# Authors:      Jeff Terpstra                      Joe McKean
# 
# Address:      Department of Statistics           Department of Statistics
#               North Dakota State University      Western Michigan University
#               P.O. Box 5575 Waldron 201H         5506 Everett Tower
#               Fargo, ND  58105                   Kalamazoo, MI  49008
#
# Email:        Jeff.Terpstra@ndsu.edu             joseph.mckean@wmich.edu
#
# Phone:        701-231-8188                       269-387-4541
#
# Fax:          701-231-8734                       269-387-1419                     

# Function:
#  wwfit
#
#Description:  
#  The function wwfit is used to calculate weighted Wilcoxon (WW) estimates.
#
#Usage: 
#  wwfit(x, y, bij=wilwts(as.matrix(x)), center=F) 
#
#Arguments:
#       x: n x p design matrix (without an initial column of ones for the
#          intercept) where the columns represent variables. 
#       y: n x 1 vector of responses
#     bij: n(n-1)/2 x 1 vector of weights to be used for the (i,j)th
#          comparison.  bij can be calculated using the functions wilwts,
#          grwts, and hbrwts; or one can supply other user defined weights.
#  center: a logical value.  If center = T then the design matrix is 
#          centered with the sample mean first.  If center = F then the original 
#          design matrix is not centered.
# 
#Value:
#  list(coefficients=c(int,est), residuals=resid, weights=wts)
#
#  coefficients: int = median-based residual estimate of the intercept parameter.
#                est = ww-estimates of the p regression parameters.
#     residuals: ww-estimate residuals.  Note, median(residuals)=0.
#       weights: n x n symmetric matrix of weights (bij).  The diagonal is set to 0. 
#
#Details:
#  wwfit minimizes a weighted version of Gini's mean difference dispersion function.
#  In particular, it can be used to calculate the WIL, GR, and HBR-estimates.
#  The algorithm is essentially a weighted L1 algorithm with O(n^2) terms.  This
#  is analogous to what is done in weighted least squares regression.  This function
#  requires the "quantreg" package which is available from http://cran.r-project.org/.
#
#References:
#  Hettmansperger, T. P. and McKean, J. W. Robust nonparametric statistical methods.
#    Edward Arnold, London, 1998.
#
#  Sievers, Gerald L. (1983). A weighted dispersion function for estimation in linear
#    models. Comm. Statist. A---Theory Methods, 12(10), 1161--1179.
#
#See Also:
#  pairup, wilwts, grwts, hbrwts
#
#Examples:
#
#> # Wilcoxon fit of telephone data from HM(Example 3.3.1, p.151).
#> calls
#   Year Number
#1    50   0.44
#2    51   0.47
#.     .      .
#24   73   2.90
#
#> wwfit(calls[,1],calls[,2])
#$coefficients
#[1] -7.1325  0.1450
#
#$residuals
# [1]  0.3225  0.2075  0.0625  0.0375 -0.0375 -0.1125 -0.1775 -0.2525 -0.2175
#[10] -0.2225 -0.2175 -0.2225 -0.2475  0.1175  9.7525 10.1075 11.7625 13.3175
#[19] 15.4725 18.3275  1.2825 -0.7625 -0.6075 -0.5525
#
#> # Hypothetical example illustrating use of weights (i.e. bij).
#> xy
#         x1    x2     y
# [1,]  1.20  0.36  3.71
# [2,]  0.65  1.23  4.04
# [3,]  0.68  1.53  5.02
# [4,]  0.17  0.21  2.66
# [5,] -0.69  0.66  1.00
# [6,]  1.18  1.26  3.65
# [7,]  0.30 -1.07 -0.17
# [8,]  0.79 -0.37  2.52
# [9,] -0.27 -0.35  0.97
#[10,]  0.56  0.36  1.46
#[11,] -1.59  0.89  1.78
#[12,]  0.59 -0.65  0.11
#[13,]  1.82  0.81  2.51
#
#> # All weights equal to 1.
#> wwfit(xy[,1:2],xy[,3])$coef
#               x1       x2 
#1.634837 0.712744 1.436433 
#
#> # All weights not equal to 1.
#> h=pairup(c(rep(.5,2),rep(1,11)))
#> wwfit(xy[,1:2],xy[,3],bij=h[,1]*h[,2])
#$coefficients
#                 x1        x2 
#1.3132927 0.5152459 1.4448857 
#
#$residuals
# [1]  1.258253e+00  6.145881e-01  1.145665e+00  9.556895e-01 -9.113976e-01 -9.183878e-02
# [7] -9.183878e-02  1.334271e+00  3.015337e-01 -6.619892e-01 -1.110223e-16 -5.681121e-01
#[13] -9.113976e-01
#
#$weights
#      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
# [1,] 0.00 0.25  0.5  0.5  0.5  0.5  0.5  0.5  0.5   0.5   0.5   0.5   0.5
# [2,] 0.25 0.00  0.5  0.5  0.5  0.5  0.5  0.5  0.5   0.5   0.5   0.5   0.5
# [3,] 0.50 0.50  0.0  1.0  1.0  1.0  1.0  1.0  1.0   1.0   1.0   1.0   1.0
# [4,] 0.50 0.50  1.0  0.0  1.0  1.0  1.0  1.0  1.0   1.0   1.0   1.0   1.0
# [5,] 0.50 0.50  1.0  1.0  0.0  1.0  1.0  1.0  1.0   1.0   1.0   1.0   1.0
# [6,] 0.50 0.50  1.0  1.0  1.0  0.0  1.0  1.0  1.0   1.0   1.0   1.0   1.0
# [7,] 0.50 0.50  1.0  1.0  1.0  1.0  0.0  1.0  1.0   1.0   1.0   1.0   1.0
# [8,] 0.50 0.50  1.0  1.0  1.0  1.0  1.0  0.0  1.0   1.0   1.0   1.0   1.0
# [9,] 0.50 0.50  1.0  1.0  1.0  1.0  1.0  1.0  0.0   1.0   1.0   1.0   1.0
#[10,] 0.50 0.50  1.0  1.0  1.0  1.0  1.0  1.0  1.0   0.0   1.0   1.0   1.0
#[11,] 0.50 0.50  1.0  1.0  1.0  1.0  1.0  1.0  1.0   1.0   0.0   1.0   1.0
#[12,] 0.50 0.50  1.0  1.0  1.0  1.0  1.0  1.0  1.0   1.0   1.0   0.0   1.0
#[13,] 0.50 0.50  1.0  1.0  1.0  1.0  1.0  1.0  1.0   1.0   1.0   1.0   0.0
# 
#Function:
#  wwest
#
#Description:
#  This function performs a Weighted Wilcoxon analysis using (the default)
#  WIL, GR, or HBR weights.  If bij is numeric then GR (i.e. non-random)
#  weights are assumed.
#
#Usage:
#  wwest(x, y, bij="WIL", center=F, print.tbl=T) 
#
#Arguments:
#          x: n x p design matrix (without an initial column of ones for the
#             intercept) where the columns represent variables. 
#          y: n x 1 vector of responses
#        bij: One of "WIL", "GR", or "HBR".  bij can also be a  
#             n(n-1)/2 x 1 vector of weights in which case wwest performs a 
#             GR analysis using bij.
#     center: a logical value.  If center = T then the design matrix is 
#             centered with the sample mean first.  If center = F then the original 
#             design matrix is not centered.
#  print.tbl: a logical value.  If print.tbl = T then the test results for
#             H0: Beta_j = 0, j=1,2,...,p, are printed along with the results
#             of a test for regression significance.  If print.tbl = F the table
#             is returned in a list.
#
#Value:
#  invisible(list(tmp1=tmp1,tmp2=tmp2,ans=ans))
#
#  tmp1: the returned list from wwfit.  
#  tmp2: the returned list from varcov.gr or varcov.hbr, which contains the 
#        estimated variance-covariance matrix for the parameter estimates.
#   ans: a hypothesis testing table corresponding to H0: Beta_j = 0, j=1,2,...,p.
#
#Details:
#  wwest essentially returns inference results for a "WIL", "GR", and or a "HBR" fit
#  of the data.  If print.tbl = T then the results of a Wald-based test of 
#  H0: Beta_1 = Beta_2 = ... = Beta_p = 0 are printed along with the results of
#  Wald-based tests of H0: Beta_j = 0, j=1,2,...,p.  A graphical residual analysis
#  is also available if print.tbl = T.  If print.tbl = F then a list of lists is
#  returned.  This master list contains information (e.g. estimates, residuals, 
#  weights, variance-covariance matrices, scale parameter estimates, test statistics,
#  p-values, etc.) from the functions wwfit, varcov.gr, varcov.hbr, and wwest itself.
#
#References:
#  Hettmansperger, T. P. and McKean, J. W. Robust nonparametric statistical methods.
#    Edward Arnold, London, 1998.
#
#  Sievers, Gerald L. (1983). A weighted dispersion function for estimation in linear
#    models. Comm. Statist. A---Theory Methods, 12(10), 1161--1179.
#
#See Also:
#  wwfit, pairup, wilwts, grwts, hbrwts, varcov.gr, varcov.hbr, studres.gr,
#  studres.hbr, wald
#
#Examples:
#
#> # GR fit of stars data of HM(Example 5.3.1, p.289).
#> wwest(stars[,1],stars[,2],"GR")
#
#Wald Test of H0: BETA1=0
#TS: 36.9845 PVAL: 0 
#
#          EST     SE    TVAL  PVAL
#BETA0 -6.9800 1.9351 -3.6071 8e-04
#BETA1  2.7273 0.4485  6.0815 0e+00
#
#Would you like to see residual plots (y/n)? 
#n
#
#> # HBR fit of baseball salaries data from HM(Example 3.3.2, p.152).
#> wwest(baseball[,1:7],baseball[,8],"HBR")
#
#Wald Test of H0: BETA1=BETA2=BETA3=BETA4=BETA5=BETA6=BETA7=0
#TS: 144.0241 PVAL: 0 
#
#          EST     SE    TVAL   PVAL
#BETA0  3.9801 0.2730 14.5781 0.0000
#BETA1  0.8678 0.0396 21.9103 0.0000
#BETA2  0.0471 0.0239  1.9714 0.0503
#BETA3 -0.0430 0.0242 -1.7754 0.0776
#BETA4 -0.1065 0.0556 -1.9156 0.0571
#BETA5  0.0080 0.0031  2.5554 0.0115
#BETA6  0.0052 0.0025  2.1296 0.0347
#BETA7  0.0096 0.0102  0.9369 0.3502
#
#Would you like to see residual plots (y/n)? 
#n
#
#Function:
#  pairup
#
#Description:
#  This function sets up pairwise comparisons.
#
#Usage:
#  pairup(x, type="less")
#
#Arguments:
#     x: a n x 1 vector or n x p matrix.
#  type: one of "less", "leq", or "neq".
#
#Value:
#  A matrix which contains the (row based) pairwise comparisons of x.  If type="less"
#  the n(n-1)/2 (i.e. i<j) pairwise comparisons are returned.  If type="leq" the 
#  n(n-1)/2 + n (i.e. i<=j) pairwise comparisons are returned.  If type="neq" the
#  n*n (i.e. all) pairwise comparisons are returned.
#
#Details:
#  The function pairup takes a n x 1 vector or n x p matrix as input and 
#  constructs a new matrix which contains the pairwise comparisons.  That is,
#  each row of the returned value corresponds to a pairwise comparison.  The 
#  number of pairwise comparisons returned is controled by the "type" argument.
#
#References:
#
#See Also:
#  rep, apply
#
#Examples:
#
#> pairup(1:4)
#     [,1] [,2]
#[1,]    1    2
#[2,]    1    3
#[3,]    1    4
#[4,]    2    3
#[5,]    2    4
#[6,]    3    4
#
#> x=cbind(c(11,21,31,41),c(12,22,32,42),c(13,23,33,43))
#> x
#     [,1] [,2] [,3]
#[1,]   11   12   13
#[2,]   21   22   23
#[3,]   31   32   33
#[4,]   41   42   43
#
#> pairup(x,type="leq")
#      [,1] [,2] [,3] [,4] [,5] [,6]
# [1,]   11   12   13   11   12   13
# [2,]   11   12   13   21   22   23
# [3,]   11   12   13   31   32   33
# [4,]   11   12   13   41   42   43
# [5,]   21   22   23   21   22   23
# [6,]   21   22   23   31   32   33
# [7,]   21   22   23   41   42   43
# [8,]   31   32   33   31   32   33
# [9,]   31   32   33   41   42   43
#[10,]   41   42   43   41   42   43 
#
#Function:
#  wilwts
#
#Description:
#  Wilcoxon weights.
#
#Usage:
#  wilwts(xmat) 
#
#Arguments:
#  xmat: a n x p design matrix.
#
#Value:
#  A (n(n-1)/2) x 1 vector of ones.
#
#Details:
#  wilwts simply returns a n(n-1)/2 x 1 vector of ones.  These weights correspond to
#  the well-known Wilcoxon estimate.
#
#References:
#  Hettmansperger, T. P. and McKean, J. W. Robust nonparametric statistical methods.
#  Edward Arnold, London, 1998.
#
#See Also:
#  
#Examples:
#
#> x
#     [,1] [,2] [,3]
#[1,]   11   12   13
#[2,]   21   22   23
#[3,]   31   32   33
#[4,]   41   42   43
#
#> wilwts(x)
#[1] 1 1 1 1 1 1
#
#Function:
#  theilwts
#
#Description:
#  Theil weights.
#
#Usage:
#  theilwts(xmat) 
#
#Arguments:
#  xmat: a n x p design matrix.
#
#Value:
#  The (n(n-1)/2) x 1 vector 1/||x_j - x_i|| where ||z|| denotes the Euclidean
#  norm of z.
#
#Details:
#  The function "theilwts" calculates and returns weights that correspond to the
#  theil type estimate.  Specifically, the following weights are returned
#
#  b_{ij} = 1/||x_j - x_i|| 
#
#  where || || represents the Euclidean norm.  If the denominator is zero then
#  b_{ij} is set to zero.  Note that in the simple regression case this reduces
#  to 1/|x_j - x_i| and the theil estimate corresponds to the median of the
#  sample slopes.
#
#References:
#  Daniel, W. W. Applied Nonparametric Statistics (Second Edition). PWS-KENT
#  Publishing Company, Wadsworth Inc., Boston, 1990.
#    
#See Also:
#  
#Examples:
#
#> x
#     [,1] [,2] [,3]
#[1,]   11   12   13
#[2,]   21   22   23
#[3,]   31   32   33
#[4,]   41   42   43
#
#> round(theilwts(x),4)
#[1] 0.0577 0.0289 0.0192 0.0577 0.0289 0.0577
#
#Function:
#  grwts
#
#Description:
#  Weights which correspond to the generalized rank (GR) estimate of 
#  Naranjo and Hettmansperger (1994) and Chang et al. (1999).
# 
#Usage:
#  grwts(xmat, robdis2=mycov.rob(as.matrix(xmat),method="mcd")$robdis2,
#        percent=0.95, k=2)
#
#Arguments:
#     xmat: a n x p design matrix.
#  robdis2: a n x 1 vector of (squared) robust-type Mahalanobis distances 
#           corresponding to the rows of xmat.  The default is based on 
#           the Minimum Covariance Determinant (MCD) estimate of location
#           and dispersion.
#  percent: denotes the percentile of a chi-square(p) distribution; used to 
#           determine cutoff values for outlying observations in the design
#           matrix.
#        k: denotes a tuning parameter corresponding to an exponent (k/2) that
#           controls the severity of downweighting.  When k=0 the Wilcoxon
#           estimate is obtained.  Higher values of k indicate severe
#           downweighting. 
#
#Value:
#  A (n(n-1)/2) x 1 vector of (factored) weights.
#
#Details:
#  The function "grwts" takes as input a n x p design matrix, a n x 1 vector of
#  (squared) robust distances, and some tuning constants and returns a
#  (n(n-1)/2) x 1 vector of weights; bij = hi * hj where hi is defined as in
#  Chang et al. (1999), pg. 207.  These are the weights that correspond to the
#  so-called GR-estimate.
#
#References:
#  Chang, W.H., McKean, J.W., Naranjo, J.D., and Sheather, S.J. (1999). "High
#  Breakdown Rank Regression". Journal of the American Statistical Association,
#  94(445), 205-219.
#
#  Naranjo, J.D. and Hettmansperger, T.P. (1994). "Bounded Influence Rank Regression".
#  Journal of the Royal Statistical Society, Series B, 56(1), 209-220.
#
#  Hettmansperger, T. P. and McKean, J. W. Robust nonparametric statistical methods.
#  Edward Arnold, London, 1998.
#
#See Also:
#  [my]cov.rob, pairup
#
#Examples:
#
#> xy
#         x1    x2     y
# [1,]  1.20  0.36  3.71
# [2,]  0.65  1.23  4.04
# [3,]  0.68  1.53  5.02
# [4,]  0.17  0.21  2.66
# [5,] -0.69  0.66  1.00
# [6,]  1.18  1.26  3.65
# [7,]  0.30 -1.07 -0.17
# [8,]  0.79 -0.37  2.52
# [9,] -0.27 -0.35  0.97
#[10,]  0.56  0.36  1.46
#[11,] -1.59  0.89  1.78
#[12,]  0.59 -0.65  0.11
#[13,]  1.82  0.81  2.51
#
#> # Based on MCD
#> grwts(xy[,1:2])
# [1] 1.00000000 1.00000000 1.00000000 0.48455158 1.00000000 1.00000000 1.00000000 1.00000000
# [9] 1.00000000 0.16893383 1.00000000 0.72303264 1.00000000 1.00000000 0.48455158 1.00000000
#[17] 1.00000000 1.00000000 1.00000000 1.00000000 0.16893383 1.00000000 0.72303264 1.00000000
#[25] 0.48455158 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 0.16893383 1.00000000
#[33] 0.72303264 0.48455158 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 0.16893383
#[41] 1.00000000 0.72303264 0.48455158 0.48455158 0.48455158 0.48455158 0.48455158 0.08185715
#[49] 0.48455158 0.35034661 1.00000000 1.00000000 1.00000000 1.00000000 0.16893383 1.00000000
#[57] 0.72303264 1.00000000 1.00000000 1.00000000 0.16893383 1.00000000 0.72303264 1.00000000
#[65] 1.00000000 0.16893383 1.00000000 0.72303264 1.00000000 0.16893383 1.00000000 0.72303264
#[73] 0.16893383 1.00000000 0.72303264 0.16893383 0.12214467 0.72303264
#
#> # Based on MVE
#> grwts(xy[,1:2],robdis2=mycov.rob(as.matrix(xy[,1:2]),method="mve")$robdis2)
# [1] 1.00000000 1.00000000 1.00000000 0.48455158 1.00000000 1.00000000 1.00000000 1.00000000
# [9] 1.00000000 0.16893383 1.00000000 0.72303264 1.00000000 1.00000000 0.48455158 1.00000000
#[17] 1.00000000 1.00000000 1.00000000 1.00000000 0.16893383 1.00000000 0.72303264 1.00000000
#[25] 0.48455158 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 0.16893383 1.00000000
#[33] 0.72303264 0.48455158 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 0.16893383
#[41] 1.00000000 0.72303264 0.48455158 0.48455158 0.48455158 0.48455158 0.48455158 0.08185715
#[49] 0.48455158 0.35034661 1.00000000 1.00000000 1.00000000 1.00000000 0.16893383 1.00000000
#[57] 0.72303264 1.00000000 1.00000000 1.00000000 0.16893383 1.00000000 0.72303264 1.00000000
#[65] 1.00000000 0.16893383 1.00000000 0.72303264 1.00000000 0.16893383 1.00000000 0.72303264
#[73] 0.16893383 1.00000000 0.72303264 0.16893383 0.12214467 0.72303264
#
#Function:
#  hbrwts
#
#Description:
#  Weights which correspond to the high breakdown rank (HBR) estimate of 
#  Chang et al. (1999).
#
#Usage:
#  hbrwts(xmat, y, robdis2=mycov.rob(as.matrix(xmat),method="mcd")$robdis2,
#         percent=0.95, intest=myltsreg(xmat,y)$coef)
#
#Arguments:
#     xmat: a n x p design matrix.
#        y: a n x 1 vector of responses.
#  robdis2: a n x 1 vector of (squared) robust-type Mahalanobis distances 
#           corresponding to the rows of xmat.  The default is based on 
#           the Minimum Covariance Determinant (MCD) estimate of location
#           and dispersion.
#  percent: denotes the percentile of a chi-square(p) distribution; used to 
#           determine cutoff values for outlying observations in the design
#           matrix.  The default is 0.95.
#   intest: a p x 1 vector of initial regression parameter estimates.  The 
#           default uses Least Trimmed (50%) Squares (LTS) estimates.  These 
#           estimates are used to calculate an initial set of residuals; which
#           determine the weights.
# 
#Value:
#  A (n(n-1)/2) x 1 vector of (non-factored) high breakdown rank weights as
#  defined in Chang et al. (1999).
#
#Details:
#  The function "hbrwts" takes as input a n x p design matrix, a n x 1 vector of
#  responses, a n x 1 vector of (squared) robust distances, a p x 1 vector of
#  initial regression parameter estimates, and a tuning constant and returns a
#  (n(n-1)/2) x 1 vector of weights, say bij, where bij is defined as in
#  Chang et al., pg. 206.  These are the weights that correspond to the so-called
#  HBR-estimate.
# 
#References:
#  Chang, W.H., McKean, J.W., Naranjo, J.D., and Sheather, S.J. (1999). "High
#  Breakdown Rank Regression". Journal of the American Statistical Association,
#  94(445), 205-219.
#
#  Hettmansperger, T. P. and McKean, J. W. Robust nonparametric statistical methods.
#  Edward Arnold, London, 1998.
#
#  Rousseeuw, P. J. and Leroy, A. M.  Robust Regression and Outlier Detection.
#  Wiley, 1987.
#
#  Rousseeuw, P. J. and van Driessen, K. (1999).  "A fast algorithm for
#     the minimum covariance determinant estimator". Technometrics, 41,
#     212-223.
#
#See Also:
#  [my]cov.rob, [my]ltsreg, [my]lmsreg, mad, psi, pairup
#
#Examples:
#
#> xy
#         x1    x2     y
# [1,]  1.20  0.36  3.71
# [2,]  0.65  1.23  4.04
# [3,]  0.68  1.53  5.02
# [4,]  0.17  0.21  2.66
# [5,] -0.69  0.66  1.00
# [6,]  1.18  1.26  3.65
# [7,]  0.30 -1.07 -0.17
# [8,]  0.79 -0.37  2.52
# [9,] -0.27 -0.35  0.97
#[10,]  0.56  0.36  1.46
#[11,] -1.59  0.89  1.78
#[12,]  0.59 -0.65  0.11
#[13,]  1.82  0.81  2.51
#
#> # Based on MCD and LTS
#> hbrwts(xy[,1:2],xy[,3])
# [1] 1.00000000 1.00000000 1.00000000 0.37126866 1.00000000 1.00000000
# [7] 1.00000000 1.00000000 1.00000000 0.12045209 1.00000000 1.00000000
#[13] 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
#[19] 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
#[25] 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
#[31] 0.81388725 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
#[37] 1.00000000 1.00000000 1.00000000 0.52009519 1.00000000 1.00000000
#[43] 1.00000000 1.00000000 0.34613891 1.00000000 0.69160494 0.03674113
#[49] 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
#[55] 0.74252746 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
#[61] 0.61425320 1.00000000 1.00000000 1.00000000 1.00000000 0.11229915
#[67] 1.00000000 1.00000000 1.00000000 0.56705257 1.00000000 1.00000000
#[73] 0.22437999 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
#
#> # Based on LS estimates
#> hbrwts(xy[,1:2],xy[,3],
#+ robdis2=mycov.rob(xy[,1:2],method="classical")$robdis2,
#+ intest=lsfit(xy[,1:2],xy[,3])$coef)
# [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#[39] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#[77] 1 1
#
#Function:
#  psi
#
#Description:
#  The function "psi" takes as input a vector and returns a vector of the same
#  length whos components are -1 if x[i] <= -1, x[i] if -1 < x[i] < 1, and
#  1 if x[i] >= 1.  It is called by the function "hbrwts".
#
#Usage:
#  psi(x)
#
#Arguments:
#  x: a n x 1 numerical vector; NA and Inf's are allowed.
# 
#Value:
#  A n x 1 vector whos components are -1 if x[i] <= -1, x[i] if -1 < x[i] < 1,
#  and 1 if x[i] >= 1.
#
#Details:
#
#References:
#
#See Also:
#
#Examples:
#
#> psi(c(NA,-Inf,-10,-1,0,1,10,Inf,NA))
#[1] NA -1 -1 -1  0  1  1  1 NA
#
#> u=round(rnorm(7),2)
#> u
#[1] -2.30  0.33  0.29 -1.70 -0.96 -0.59 -1.87
#> psi(u)
#[1] -1.00  0.33  0.29 -1.00 -0.96 -0.59 -1.00
#
#Function:
#  blwts
#
#Description:
#  Factored HBR-type weights that only downweight bad leverage points.
#  These weights were used in a simulation study for autoregressive time
#  series models by Terpstra, McKean and Naranjo (2001, p.413).
#
#Usage:
#  blwts(xmat, y, robdis2=mycov.rob(as.matrix(xmat),method="mcd")$robdis2,
#        percent=0.95, k=2, intest=myltsreg(xmat,y)$coef)
#
#Arguments:
#     xmat: a n x p design matrix.
#        y: a n x 1 vector of responses.
#  robdis2: a n x 1 vector of (squared) robust-type Mahalanobis distances 
#           corresponding to the rows of xmat.  The default is based on 
#           the Minimum Covariance Determinant (MCD) estimate of location
#           and dispersion.
#  percent: denotes the percentile of a chi-square(p) distribution; used to 
#           determine cutoff values for outlying observations in the design
#           matrix.  The default is 0.95.
#        k: denotes a tuning parameter corresponding to an exponent (k/2) that
#           controls the severity of downweighting.  When k=0 the Wilcoxon
#           estimate is obtained.  Higher values of k indicate severe
#           downweighting. 
#   intest: a p x 1 vector of initial regression parameter estimates.  The 
#           default uses Least Trimmed (50%) Squares (LTS) estimates.  These 
#           estimates are used to calculate an initial set of residuals; which
#           determine the weights.
# 
#Value:
#  A (n(n-1)/2) x 1 vector of (factored) high breakdown rank weights as
#  defined in Terpstra et al. (2001).
#
#Details:
#  The function "blwts" takes as input a n x p design matrix, a n x 1 vector of
#  responses, a n x 1 vector of (squared) robust distances, a p x 1 vector of
#  initial regression parameter estimates, and two tuning constants and returns
#  a (n(n-1)/2) x 1 vector of weights, say bij, where bij = hi*hj and hi is
#  setup so that only those points with a large residual and a large distance
#  (i.e. bad leverage points) are downweighted.  All other possibilities receive
#  a weight of 1.  This corresponds to the HBR2 estimate found in Terpstra et al.
#  (2001).  
# 
#References:  
#  Terpstra, Jeffrey T.; McKean, Joseph W.; Naranjo, Joshua D. (2001). "Weighted
#  Wilcoxon estimates for autoregression".  Aust. N. Z. J. Stat., 43(4),
#  399-419.
#
#  Hettmansperger, T. P. and McKean, J. W. Robust nonparametric statistical methods.
#  Edward Arnold, London, 1998.
#
#  Rousseeuw, P. J. and Leroy, A. M.  Robust Regression and Outlier Detection.
#  Wiley, 1987.
#
#See Also:
#  hbrwts, grwts, [my]cov.rob, [my]ltsreg, [my]lmsreg
#
#Examples:
#> xy
#         x1    x2     y
# [1,]  1.20  0.36  3.71
# [2,]  0.65  1.23  4.04
# [3,]  0.68  1.53  5.02
# [4,]  0.17  0.21  2.66
# [5,] -0.69  0.66  1.00
# [6,]  1.18  1.26  3.65
# [7,]  0.30 -1.07 -0.17
# [8,]  0.79 -0.37  2.52
# [9,] -0.27 -0.35  0.97
#[10,]  0.56  0.36  1.46
#[11,] -1.59  0.89  1.78
#[12,]  0.59 -0.65  0.11
#[13,]  1.82  0.81  2.51
# 
#> # Based on MCD and LTS
#> blwts(xy[,1:2],xy[,3])
# [1] 1.00000000 1.00000000 1.00000000 0.48455158 1.00000000 1.00000000
# [7] 1.00000000 1.00000000 1.00000000 0.16893383 1.00000000 1.00000000
#[13] 1.00000000 1.00000000 0.48455158 1.00000000 1.00000000 1.00000000
#[19] 1.00000000 1.00000000 0.16893383 1.00000000 1.00000000 1.00000000
#[25] 0.48455158 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
#[31] 0.16893383 1.00000000 1.00000000 0.48455158 1.00000000 1.00000000
#[37] 1.00000000 1.00000000 1.00000000 0.16893383 1.00000000 1.00000000
#[43] 0.48455158 0.48455158 0.48455158 0.48455158 0.48455158 0.08185715
#[49] 0.48455158 0.48455158 1.00000000 1.00000000 1.00000000 1.00000000
#[55] 0.16893383 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
#[61] 0.16893383 1.00000000 1.00000000 1.00000000 1.00000000 0.16893383
#[67] 1.00000000 1.00000000 1.00000000 0.16893383 1.00000000 1.00000000
#[73] 0.16893383 1.00000000 1.00000000 0.16893383 0.16893383 1.00000000
#
#Function:
#  wts
#
#Description:
#  Wrapper function that can be used to get WIL, GR, HBR, and BL weights.
#
#Usage:
#  wts(xmat, y, type="WIL", percent=0.95, k=2,
#      robdis2=if(type!="WIL") mycov.rob(as.matrix(xmat),method="mcd")$robdis2 else NULL,
#      intest=if(type=="HBR" | type=="BL") myltsreg(xmat,y)$coef else NULL)
#
#Arguments:
#     xmat: a n x p design matrix.
#        y: a n x 1 vector of responses.
#     type: one of "WIL", "GR", "HBR", or "BL".
#  percent: denotes the percentile of a chi-square(p) distribution; used to 
#           determine cutoff values for outlying observations in the design
#           matrix.  The default is 0.95.
#        k: denotes a tuning parameter corresponding to an exponent (k/2) that
#           controls the severity of downweighting.  When k=0 the Wilcoxon
#           estimate is obtained.  Higher values of k indicate severe
#           downweighting. 
#  robdis2: a n x 1 vector of (squared) robust-type Mahalanobis distances 
#           corresponding to the rows of xmat.  The default is based on 
#           the Minimum Covariance Determinant (MCD) estimate of location
#           and dispersion.
#   intest: a p x 1 vector of initial regression parameter estimates.  The 
#           default uses Least Trimmed (50%) Squares (LTS) estimates.  These 
#           estimates are used to calculate an initial set of residuals; which
#           determine the weights.
#
#Value:
#  A (n(n-1)/2) x 1 vector of weights obtained from either wilwts, grwts, hbrwts,
#  or blwts.
#
#Details:
#  This function simply calls one of wilwts, grwts, hbrwts, or blwts.
#
#References:
#  Chang, W.H., McKean, J.W., Naranjo, J.D., and Sheather, S.J. (1999). "High
#  Breakdown Rank Regression". Journal of the American Statistical Association,
#  94(445), 205-219.
#  
#  Hettmansperger, T. P. and McKean, J. W. Robust nonparametric statistical methods.
#  Edward Arnold, London, 1998.
#
#  Naranjo, J.D. and Hettmansperger, T.P. (1994). "Bounded Influence Rank Regression".
#  Journal of the Royal Statistical Society, Series B, 56(1), 209-220.
#
#  Rousseeuw, P. J. and Leroy, A. M.  Robust Regression and Outlier Detection.
#  Wiley, 1987.
#
#  Rousseeuw, P. J. and van Driessen, K. (1999).  "A fast algorithm for
#  the minimum covariance determinant estimator". Technometrics, 41,
#  212-223.
#
#  Terpstra, Jeffrey T.; McKean, Joseph W.; Naranjo, Joshua D. (2001). "Weighted
#  Wilcoxon estimates for autoregression".  Aust. N. Z. J. Stat., 43(4),
#  399-419.
#
#See Also:
#  wilwts, grwts, hbrwts, blwts
#
#Examples:
#> xy
#         x1    x2     y
# [1,]  1.20  0.36  3.71
# [2,]  0.65  1.23  4.04
# [3,]  0.68  1.53  5.02
# [4,]  0.17  0.21  2.66
# [5,] -0.69  0.66  1.00
# [6,]  1.18  1.26  3.65
# [7,]  0.30 -1.07 -0.17
# [8,]  0.79 -0.37  2.52
# [9,] -0.27 -0.35  0.97
#[10,]  0.56  0.36  1.46
#[11,] -1.59  0.89  1.78
#[12,]  0.59 -0.65  0.11
#[13,]  1.82  0.81  2.51
#
#> round(wts(xy[,1:2],xy[,3],"WIL"),4)
# [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#[39] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#[77] 1 1
#
#> round(wts(xy[,1:2],xy[,3],"GR"),4)
# [1] 1.0000 1.0000 1.0000 0.4846 1.0000 1.0000 1.0000 1.0000 1.0000 0.1689
#[11] 1.0000 0.7230 1.0000 1.0000 0.4846 1.0000 1.0000 1.0000 1.0000 1.0000
#[21] 0.1689 1.0000 0.7230 1.0000 0.4846 1.0000 1.0000 1.0000 1.0000 1.0000
#[31] 0.1689 1.0000 0.7230 0.4846 1.0000 1.0000 1.0000 1.0000 1.0000 0.1689
#[41] 1.0000 0.7230 0.4846 0.4846 0.4846 0.4846 0.4846 0.0819 0.4846 0.3503
#[51] 1.0000 1.0000 1.0000 1.0000 0.1689 1.0000 0.7230 1.0000 1.0000 1.0000
#[61] 0.1689 1.0000 0.7230 1.0000 1.0000 0.1689 1.0000 0.7230 1.0000 0.1689
#[71] 1.0000 0.7230 0.1689 1.0000 0.7230 0.1689 0.1221 0.7230
#
#> round(wts(xy[,1:2],xy[,3],"HBR"),4)
# [1] 1.0000 1.0000 1.0000 0.3713 1.0000 1.0000 1.0000 1.0000 1.0000 0.1205
#[11] 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000
#[21] 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000
#[31] 0.8139 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.5201
#[41] 1.0000 1.0000 1.0000 1.0000 0.3461 1.0000 0.6916 0.0367 1.0000 1.0000
#[51] 1.0000 1.0000 1.0000 1.0000 0.7425 1.0000 1.0000 1.0000 1.0000 1.0000
#[61] 0.6143 1.0000 1.0000 1.0000 1.0000 0.1123 1.0000 1.0000 1.0000 0.5671
#[71] 1.0000 1.0000 0.2244 1.0000 1.0000 1.0000 1.0000 1.0000
#
#> round(wts(xy[,1:2],xy[,3],"BL"),4)
# [1] 1.0000 1.0000 1.0000 0.4846 1.0000 1.0000 1.0000 1.0000 1.0000 0.1689
#[11] 1.0000 1.0000 1.0000 1.0000 0.4846 1.0000 1.0000 1.0000 1.0000 1.0000
#[21] 0.1689 1.0000 1.0000 1.0000 0.4846 1.0000 1.0000 1.0000 1.0000 1.0000
#[31] 0.1689 1.0000 1.0000 0.4846 1.0000 1.0000 1.0000 1.0000 1.0000 0.1689
#[41] 1.0000 1.0000 0.4846 0.4846 0.4846 0.4846 0.4846 0.0819 0.4846 0.4846
#[51] 1.0000 1.0000 1.0000 1.0000 0.1689 1.0000 1.0000 1.0000 1.0000 1.0000
#[61] 0.1689 1.0000 1.0000 1.0000 1.0000 0.1689 1.0000 1.0000 1.0000 0.1689
#[71] 1.0000 1.0000 0.1689 1.0000 1.0000 0.1689 0.1689 1.0000
#
#Function:
#  mymahalanobis
#
#Description:
#  Modified version of mahalanobis.
#
#Usage:
#  mymahalanobis(x, center, cov, inverted = FALSE, tol.inv = 1e-07)
#
#Arguments:
#       x: vector or matrix of data with, say, p columns.
#  center: mean vector of the distribution or second data vector of
#          length p.
#     cov: covariance matrix (p x p) of the distribution.
#inverted: logical.  If `TRUE', `cov' is supposed to contain the inverse
#          of the covariance matrix.
# tol.inv: tolerance to be used for computing the inverse (if `inverted'
#          is false), see `solve'.
#
#Value:
#  Returns the Mahalanobis distance of all rows in x and the vector
#  mu=center with respect to Sigma=cov. This is (for vector x)
#  defined as: D^2 = (x - mu)' Sigma^{-1} (x - mu).
#
#Details:
#  This function is identical to mahalanobis except that the matrix statement
#  in the third row has been changed from matrix(x, ncol = length(x)) to
#  matrix(x, nrow = length(x)) so that the function works for p=1. 
#
#References:
#
#See Also:
#  mahalanobis
#
#Examples:
#> x=round(rnorm(10),2)
#> x
# [1]  0.15  0.01  0.48  1.33 -2.10  1.08 -1.23 -0.18  0.45 -1.31
#
#> mahalanobis(x,mean(x),var(x))
#Error in x %*% cov : non-conformable arguments
#
#> mymahalanobis(x,mean(x),var(x))
# [1] 0.065734628 0.016667585 0.309598492 1.766813618 3.201446001 1.214230765
# [7] 0.996553624 0.001904489 0.279989640 1.147061157
#
#Function:
#  mycov.rob
#
#Description:
#  Modified version of cov.rob.
#
#Usage:
#  mycov.rob(x, cor = FALSE, quantile.used = floor((n + p + 1)/2), 
#            method = c("mve", "mcd", "classical"), nsamp = "best")
#
#Arguments:
#  See help(cov.rob).
#
#Value:
#  See help(cov.rob).
#
#Details:
#  This function is identical to cov.rob except for the following changes.
#  1. This version eliminates the use of "seed" and sets the seed to the same
#     value every time. This is similar to what S+ does.
#  2. This version uses mymahalanobis instead of mahalanobis.
#  3. This version returns squared (robust) distances, robdis2.
#
#References:
#  See help(cov.rob).
#
#See Also:
#  cov.rob
#
#Examples:
#> data(stackloss)
#
#> cov.rob(stackloss)
#$center
#  Air.Flow Water.Temp Acid.Conc. stack.loss 
#   56.3750    20.0000    85.4375    13.0625 
#
#$cov
#            Air.Flow Water.Temp Acid.Conc. stack.loss
#Air.Flow   23.050000   6.666667  16.625000  19.308333
#Water.Temp  6.666667   5.733333   5.333333   7.733333
#Acid.Conc. 16.625000   5.333333  34.395833  13.837500
#stack.loss 19.308333   7.733333  13.837500  18.462500
#
#$msg
#[1] "27 singular samples of size 5 out of 2500"
#
#$crit
#[1] 19.39599
#
#$best
# [1]  5  6  7  8  9 13 14 15 16 17 18 19 20
#
#$n.obs
#[1] 21
#
#> mycov.rob(stackloss)
#$center
#  Air.Flow Water.Temp Acid.Conc. stack.loss 
#   56.3750    20.0000    85.4375    13.0625 
#
#$cov
#            Air.Flow Water.Temp Acid.Conc. stack.loss
#Air.Flow   23.050000   6.666667  16.625000  19.308333
#Water.Temp  6.666667   5.733333   5.333333   7.733333
#Acid.Conc. 16.625000   5.333333  34.395833  13.837500
#stack.loss 19.308333   7.733333  13.837500  18.462500
#
#$robdis2
#         1          2          3          4          5          6          7 
#105.206242  47.174991  84.653702  76.305096   1.700333   2.523399   3.308281 
#         8          9         10         11         12         13         14 
#  4.132534   2.475907   4.174699   3.022447   3.737464   6.671697   4.678897 
#        15         16         17         18         19         20         21 
#  5.007879   3.249796   5.658685   2.172216   2.754463   4.731303  62.836886 
#
#$msg
#[1] "26 singular samples of size 5 out of 2500"
#
#$crit
#[1] 19.39599
#
#$best
# [1]  5  6  7  8  9 13 14 15 16 17 18 19 20
#
#$n.obs
#[1] 21
#
#Function:
#  myltsreg
#
#Description:
#  Modified version of ltsreg.
#
#Usage:
#  myltsreg(xmat, y)
# 
#Arguments:
#  xmat: a n x p matrix or data frame containing the explanatory variables,
#        without an intercept column.
#     y: a n x 1 response vector.
#
#Value:
#  list(coefficients=tmp$coefficients, residuals=tmp$residuals)
#
#  coefficients: (p+1) x 1 vector of regression parameter estimates.
#     residuals: n x 1 vector of residuals.
#
#Details:
#  This function is essentially the same as ltsreg except for the following changes.
#  1. This version eliminates the use of "seed" and sets the seed to the same
#     value every time. This is similar to what S+ does.
#  2. This version only returns the coefficients and the residuals.
#  3. This version only uses the default settings for the various options in 
#     ltsreg.
#
#References:
#  See help(ltsreg).
#
#See Also:
#  ltsreg
#
#Examples:
#> data(stackloss)
#
#> ltsreg(stackloss[,1:3],stackloss[,4],intercept=T)$coef
#  (Intercept)      Air.Flow    Water.Temp    Acid.Conc. 
#-3.580556e+01  7.500000e-01  3.333333e-01 -2.105723e-17 
#
#> myltsreg(stackloss[,1:3],stackloss[,4])$coef
#  (Intercept)      Air.Flow    Water.Temp    Acid.Conc. 
#-3.580556e+01  7.500000e-01  3.333333e-01 -1.183191e-17 
#
#Function:
#  mylmsreg
#
#Description:
#  Modified version of lmsreg.
#
#Usage:
#  mylmsreg(xmat, y) 
#
#Arguments:
#  xmat: a n x p matrix or data frame containing the explanatory variables,
#        without an intercept column.
#     y: a n x 1 response vector.
#
#Value:
#  list(coefficients=tmp$coefficients, residuals=tmp$residuals)
#
#  coefficients: (p+1) x 1 vector of regression parameter estimates.
#     residuals: n x 1 vector of residuals.
#
#Details:
#  This function is essentially the same as lmsreg except for the following changes.
#  1. This version eliminates the use of "seed" and sets the seed to the same
#     value every time. This is similar to what S+ does.
#  2. This version only returns the coefficients and the residuals.
#  3. This version only uses the default settings for the various options in 
#     lmsreg.
#
#References:
#  See help(lmsreg).
#
#See Also:
#  lmsreg
#
#Examples:
#> data(stackloss)
#
#> lmsreg(stackloss[,1:3],stackloss[,4],intercept=T)$coef
#  (Intercept)      Air.Flow    Water.Temp    Acid.Conc. 
#-3.425000e+01  7.142857e-01  3.571429e-01 -2.351731e-16 
#
#> mylmsreg(stackloss[,1:3],stackloss[,4])$coef
#  (Intercept)      Air.Flow    Water.Temp    Acid.Conc. 
#-3.333333e+01  7.083333e-01  3.333333e-01  1.477219e-16 
#
#> lmsreg(stackloss[,1:3],stackloss[,4],intercept=T)$coef
#  (Intercept)      Air.Flow    Water.Temp    Acid.Conc. 
#-3.550000e+01  5.000000e-01  1.000000e+00  1.872808e-17 
#
#> mylmsreg(stackloss[,1:3],stackloss[,4])$coef
#  (Intercept)      Air.Flow    Water.Temp    Acid.Conc. 
#-3.333333e+01  7.083333e-01  3.333333e-01  1.477219e-16 
#
#> lmsreg(stackloss[,1:3],stackloss[,4],intercept=T)$coef
#  (Intercept)      Air.Flow    Water.Temp    Acid.Conc. 
#-3.425000e+01  7.142857e-01  3.571429e-01 -7.658801e-17 
#
#> mylmsreg(stackloss[,1:3],stackloss[,4])$coef
#  (Intercept)      Air.Flow    Water.Temp    Acid.Conc. 
#-3.333333e+01  7.083333e-01  3.333333e-01  1.477219e-16 
#
#> # Note that lmsreg can change from call to call, but mylmsreg stays 
#> # the same.
#
#Function:
#  wilcoxontau
#
#Description:
#  Density estimate of 1/(sqrt(12)*\int f^2(x)dx) discussed in
#  Section 3.7.1 of Hettmansperger and McKean (1998).
#
#Usage:
#  wilcoxontau(resd,p,delta=if((length(resd)/p) > 5) 0.80 else 0.95,param=2)
#
#Arguments:
#  resd: a n x 1 vector of (full model) residuals.
#     p: the number of regression coefficients (without the intercept).
# delta: bandwidth parameter for estimate of wilcoxontau.  Values should be 
#        between .8 and .99.  Larger values result in a more conservative test.
#        If n/p is less than 5, recommend delta = .95.  Default is .8.
# param: a parameter for Huber's degrees of freedom correcton
#        as discussed in Huber (1981, p.174).  It is used in the estimator
#        of wilcoxontau.  Smaller values lead to more conservative estimates
#        of wilcoxontau.  The default is 2.  This is the benchmark for
#        labeling a standardized residual as a potential outlier.
#
#Value:
#  A scaler which estimates 1/(sqrt(12)*\int f^2(x)dx).
# 
#Details:
#  See Section 3.7.1 of Hettmansperger and Mckean (1998) and or 
#  Koul, Sievers, and McKean (1987).
#
#References:
#  Hettmansperger, T. P. and McKean, J. W. Robust nonparametric statistical methods.
#  Edward Arnold, London, 1998.
#
#  Huber, Peter J. Robust statistics. Wiley Series in Probability and
#  Mathematical Statistics. John Wiley & Sons, Inc., New York, 1981.
#
#  Koul, Hira L.; Sievers, Gerald L. and McKean, Joseph. (1987). An estimator of the
#  scale parameter for the rank analysis of linear models under general score
#  functions. Scand. J. Statist. 14(2), 131-141.
#
#See Also:
#  taustar, varcov.gr, varcov.hbr, pairup
#
#Examples:
#> w=wwfit(stars[,1],stars[,2])
#
#> w$resid
# [1]  1.100000e-01  7.105607e-01 -2.424299e-01  7.105607e-01  3.663551e-02
# [6]  3.828972e-01 -7.226168e-01  2.453271e-01  3.975701e-01 -4.440892e-16
#[11]  1.905607e-01  3.585981e-01  3.524299e-01 -1.241589e+00 -8.981308e-01
#[16] -5.161682e-01 -1.246729e+00 -9.161682e-01 -1.006729e+00  3.505607e-01
#[21] -7.781308e-01 -9.381308e-01 -6.761682e-01 -2.128037e-01 -9.523364e-02
#[26] -4.361682e-01 -4.981308e-01 -2.152336e-01 -8.014953e-01  5.057944e-01
#[31] -6.952336e-01  7.056075e-02  1.381308e-01  7.505607e-01 -8.467290e-01
#[36]  6.191589e-01  5.626168e-02  1.381308e-01  1.362617e-01  4.785981e-01
#[41] -4.952336e-01 -2.186916e-02  2.819626e-01  2.581308e-01  5.057944e-01
#[46] -1.018692e-01 -5.961682e-01
#
#> wilcoxontau(w$resid,1)
#      80% 
#0.6043436 
#
#Function:
#  taustar
#
#Description:
#  Confidence interval estimate of taustar = 1/(2*f(0)) discussed on
#  pages 7-8 and 25-26 of Hettmansperger and McKean (1998).
#
#Usage:
#  taustar(resid,p,conf=.95)
#
#Arguments:
#  resid: a n x 1 vector of (full model) residuals.
#      p: the number of regression coefficients (without the intercept).
#   conf: confidence coefficient for the nonparametric confidence interval
#         being used to calculate the scale parameter estimate.
#         Default is 0.95.
#
#Value:
#  A scaler which estimates 1/(2*f(0)) where f denotes the pdf of the 
#  errors. 
# 
#Details:
#  See pages 7-8 and 25-26 of Hettmansperger and Mckean (1998) and or 
#  McKean and Schrader (1984).
#
#References:
#  Hettmansperger, T. P. and McKean, J. W. Robust nonparametric statistical methods.
#  Edward Arnold, London, 1998.
#
#  McKean, Joseph W. and Schrader, Ronald M. (1984). A comparison of methods for
#  Studentizing the sample median. Communications in Statistics,
#  Part B -- Simulation and Computation, 13, 751-773.
# 
#See Also:
#  wilcoxontau, varcov.gr, varcov.hbr
#
#Examples:
#> w=wwfit(stars[,1],stars[,2])
#
#> w$resid
# [1]  1.100000e-01  7.105607e-01 -2.424299e-01  7.105607e-01  3.663551e-02
# [6]  3.828972e-01 -7.226168e-01  2.453271e-01  3.975701e-01 -4.440892e-16
#[11]  1.905607e-01  3.585981e-01  3.524299e-01 -1.241589e+00 -8.981308e-01
#[16] -5.161682e-01 -1.246729e+00 -9.161682e-01 -1.006729e+00  3.505607e-01
#[21] -7.781308e-01 -9.381308e-01 -6.761682e-01 -2.128037e-01 -9.523364e-02
#[26] -4.361682e-01 -4.981308e-01 -2.152336e-01 -8.014953e-01  5.057944e-01
#[31] -6.952336e-01  7.056075e-02  1.381308e-01  7.505607e-01 -8.467290e-01
#[36]  6.191589e-01  5.626168e-02  1.381308e-01  1.362617e-01  4.785981e-01
#[41] -4.952336e-01 -2.186916e-02  2.819626e-01  2.581308e-01  5.057944e-01
#[46] -1.018692e-01 -5.961682e-01
#
#> taustar(w$resid,1)
#[1] 1.026483
#
#Function:
#  varcov.gr
#
#Description:
#  Variance-covariance matrix for the GR-estimate as given on page 288 of
#  Hettmansperger and McKean (1998).
#
#Usage:
#  varcov.gr(x, bmat, res, delta=0.80)
#
#Arguments:
#       x: a n x p design matrix.
#    bmat: a n x n matrix that contains the weights, bij.  This matrix is
#          returned as the "weights" list component in the function "wwfit".
#     res: a n x 1 vector of residuals.  This vector is returned as the list
#          component "residuals" in the function "wwfit".
#   delta: represents the quantile used in the estimate of wilcoxontau.  The
#          default is 0.80.
#
#Value:
#  list(varcov=varcov, tau1=tau1, tau=tau, wmat=w, cmat=cmat, vmat=vmat)
#
#  varcov: a (p+1) x (p+1) variance-covariance matrix for the GR-estimate
#          (which includes the intercept).  It corresponds to the uncentered
#          version of the design matrix.
#    tau1: a scalar which estimates 1/(2*f(0)).
#     tau: a scalar which estimates 1/(sqrt(12)*\int f^2(x)dx).
#    wmat: the weight matrix defined in (5.2.1) on page 280 of HM.  
#    cmat: the C matrix defined in (5.2.2) on page 280 of HM.  
#    vmat: the V matrix defined in (5.2.3) on page 280 of HM.  
#
#Details:
#  The function "varcov.gr" returns the variance-covariance matrix for the 
#  GR-estimate as defined on page 288 of Hettmansperger and McKean (1998).  The
#  W, C, and V matrices are those defined on page 280 of HM.  This function can
#  be used for both GR and Nonrandom (e.g. Wilcoxon) weights.
#
#References:
#  Hettmansperger, T. P. and McKean, J. W. Robust nonparametric statistical methods.
#  Edward Arnold, London, 1998.
#
#See Also:
#  wilcoxontau, taustar, varcov.hbr
#
#Examples:
#> w=wwfit(stars[,1],stars[,2],grwts(stars[,1]))
#> varcov.gr(stars[,1],w$weights,w$resid)$varcov
#          [,1]       [,2]
#[1,]  4.658024 -1.0782818
#[2,] -1.078282  0.2501814
#> sqrt(diag(varcov.gr(stars[,1],w$weights,w$resid)$varcov))
#[1] 2.1582456 0.5001813
#> # These match up with Table 5.3.2 on page 290 of HM.
#
#> w=wwfit(stars[,1],stars[,2],wilwts(stars[,1]))
#> varcov.gr(stars[,1],w$weights,w$resid)$varcov
#           [,1]       [,2]
#[1,]  1.7662552 -0.4046025
#[2,] -0.4046025  0.0938753
#> sqrt(diag(varcov.gr(stars[,1],w$weights,w$resid)$varcov))
#[1] 1.3290054 0.3063907
#> # These match up with Table 5.3.2 on page 290 of HM.
#
#Function:
#  varcov.hbr
#
#Description:
#  Variance-covariance matrix for the HBR-estimate as given in Section
#  5.5.1 of Chang et al. (1999).
#
#Usage:
#  varcov.hbr(x, bmat, res, delta=0.80)
#
#Arguments:
#       x: a n x p design matrix.
#    bmat: a n x n matrix that contains the weights, bij.  This matrix is
#          returned as the "weights" list component in the function "wwfit".
#     res: a n x 1 vector of residuals.  This vector is returned as the list
#          component "residuals" in the function "wwfit".
#   delta: represents the quantile used in the estimate of wilcoxontau.  The
#          default is 0.80.
#
#Value:
#  list(varcov=varcov, tau1=tau1, tau=tau, wmat=w, cmat=cmat, vmat=vmat)
#  
#  varcov: a (p+1) x (p+1) variance-covariance matrix for the HBR-estimate
#          (which includes the intercept).  It corresponds to the uncentered
#          version of the design matrix.
#    tau1: a scalar which estimates 1/(2*f(0)).
#     tau: a scalar which estimates 1/(sqrt(12)*\int f^2(x)dx).
#    wmat: the A matrix defined in (12) on page 210 of Chang et al.  
#    cmat: the C matrix defined in N1 on page 209 of Chang et al.  
#    vmat: the Sigma matrix defined in (10) on page 209 of Chang et al.  
#
#Details:
#  The function "varcov.hbr" returns the variance-covariance matrix for the 
#  HBR-estimate.  The W, C, and V matrices are those defined in (12), N1, and
#  (10) of Chang et al. respectively.  It corresponds to the uncentered version
#  of the design matrix.  That is, the off-diagonal elements are non-zero.
#
#References:
#  Chang, W.H., McKean, J.W., Naranjo, J.D., and Sheather, S.J. (1999). "High
#  Breakdown Rank Regression". Journal of the American Statistical Association,
#  94(445), 205-219.
#
#  Hettmansperger, T. P. and McKean, J. W. Robust nonparametric statistical methods.
#  Edward Arnold, London, 1998.
#
#See Also:
#  wilcoxontau, taustar, varcov.gr
#
#Examples:
#> # HBR
#> w=wwfit(stars[,1],stars[,2],hbrwts(stars[,1],stars[,2]))
#> varcov.hbr(stars[,1],w$weights,w$resid)$varcov
#           [,1]       [,2]
#[1,]  2.7842080 -0.6434541
#[2,] -0.6434541  0.1492933
#> sqrt(diag(varcov.hbr(stars[,1],w$weights,w$resid)$varcov))
#[1] 1.6685946 0.3863849
#
#> # GR
#> w=wwfit(stars[,1],stars[,2],grwts(stars[,1]))
#> varcov.gr(stars[,1],w$weights,w$resid)$varcov
#          [,1]       [,2]
#[1,]  4.658024 -1.0782818
#[2,] -1.078282  0.2501814
#> sqrt(diag(varcov.gr(stars[,1],w$weights,w$resid)$varcov))
#[1] 2.1582456 0.5001813
#
#> # WIL
#> w=wwfit(stars[,1],stars[,2],wilwts(stars[,1]))
#> varcov.gr(stars[,1],w$weights,w$resid)$varcov
#           [,1]       [,2]
#[1,]  1.7662552 -0.4046025
#[2,] -0.4046025  0.0938753
#> sqrt(diag(varcov.gr(stars[,1],w$weights,w$resid)$varcov))
#[1] 1.3290054 0.3063907
#
#Function:
#  wald
#
#Description:
#  This function calculates a Wald statistic and corresponding
#  p-value for the null hypothesis, H0: Ab=c, where A is a known
#  q x (p+1) matrix, b is a (p+1) x 1 vector of paramters (which
#  includes the intercept term), and c is a known (p+1) x 1 
#  vector (typically c is the zero vector).
#
#Usage:
#  wald(est, varcov, amat, true, n)
#
#Arguments:
#     est: a (p+1) x 1 estimated parameter vector which includes an
#          intercept term.  It is returned as the wwfit list component,
#          "coefficients".
#  varcov: a (p+1) x (p+1) estimated variance-covariance matrix for the
#          estimate given in est above.  It is returned as the varcov.gr
#          (or varcov.hbr) list component, "varcov".
#    amat: a q x (p+1) known matrix with rank q.
#    true: a q x 1 known vector (typically the zero vector).
#       n: the sample size.
#
#Value:
#  c(T2,pvalue)
#
#      T2: an (adjusted) Wald statistic,
#          t(Ab-c)%*%solve(A%*%Sigma%*%t(A))%*%(Ab-c)/q
#          where Sigma is the variance-covariance matrix given by varcov.
#  pvalue: P[F(q,n-p-1)>T2] where F(q,n-p-1) denotes a random variable
#          with an F-distribution with q numerator and n-p-1 denominator
#          degrees of freedom.
#
#Details:
#  This function requires that est contains both the intercept paramter 
#  estimate and the regression paramter estimates.  The use of the F
#  distribtion for calculating pvalues is documented in Section 3.6 of
#  Hettmansperger and McKean (1983, 1998) and McKean and Sheather (1991). 
#
#References:
#  Hettmansperger, T. P. and McKean, J. W. Robust nonparametric statistical methods.
#  Edward Arnold, London, 1998.
#
#  Hettmansperger, Thomas P.; McKean, Joseph W. (1983). A geometric interpretation
#  of inferences based on ranks in the linear model.  J. Amer. Statist. Assoc.
#  78 no. 384, 885-893.
#
#  McKean, Joseph W.; Sheather, Simon J. (1991). Small sample properties of robust
#  analyses of linear models based on R-estimates: a survey.  Directions in robust
#  statistics and diagnostics, Part II, New York: Springer-Verlag, 1-19.
#
#See Also:
#  varcov.gr, varcov.hbr, wwfit, wwest
#
#Examples:
#> w=wwest(baseball[,1:7],baseball[,8],"WIL")
#
#Wald Test of H0: BETA1=BETA2=BETA3=BETA4=BETA5=BETA6=BETA7=0
#TS: 114.7144 PVAL: 0 
#
#Drop Test of H0: BETA1=BETA2=BETA3=BETA4=BETA5=BETA6=BETA7=0
#TS: 57.166 PVAL: 0 
#
#          EST     SE    TVAL   PVAL
#BETA0  4.2189 0.3270 12.9002 0.0000
#BETA1  0.8390 0.0442 18.9952 0.0000
#BETA2  0.0450 0.0279  1.6142 0.1084
#BETA3 -0.0242 0.0265 -0.9138 0.3621
#BETA4 -0.1459 0.0696 -2.0969 0.0375
#BETA5  0.0061 0.0039  1.5840 0.1151
#BETA6  0.0042 0.0026  1.6099 0.1093
#BETA7  0.0121 0.0114  1.0589 0.2911
#
#> w=wwfit(baseball[,1:7],baseball[,8],wilwts(baseball[,1:7]))
#> a
#     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#[1,]    0    0    1    0    0    0    0    0
#[2,]    0    0    0    1    0    0    0    0
#[3,]    0    0    0    0    0    1    0    0
#[4,]    0    0    0    0    0    0    1    0
#[5,]    0    0    0    0    0    0    0    1
#
#> wald(w$coef,varcov.gr(baseball[,1:7],w$weights,w$resid)$varcov,a,rep(0,5),
#+ n=dim(w$weights)[1])
#[1] 1.894745e+01 6.439294e-15
#
#Function:
#  studres.gr
#
#Description:
#  This function calculates the (internally) studentized residuals
#  for the GR estimate as given on page 300 of Hettmansperger and
#  McKean (1998). See also, Naranjo et al. (1994).  These residuals
#  can be used to flag potential outliers (e.g. if the absolute 
#  value of the studentized residual exceeds 2).
#
#Usage:
#  studres.gr(x, bmat, res, delta=0.80, center=T)
#
#Arguments:
#       x: a n x p design matrix.
#    bmat: a n x n matrix that contains the weights, bij.  This matrix is
#          returned as the "weights" list component in the function "wwfit".
#     res: a n x 1 vector of residuals.  This vector is returned as the list
#          component "residuals" in the function "wwfit".
#   delta: represents the quantile used in the estimate of wilcoxontau.  The
#          default is 0.80.
#  center: should the design matrix be centered first.  The default is TRUE.
#  
#Value:
#  An n x 1 vector of studentized residuals defined as res/sqrt(diag(v)) where
#  v represents a first-order approximation of the variance-covariance matrix
#  of the raw (GR) residuals.
#
#Details:
#  The variance-covariance residual matrix estimate is based on a first-order
#  approximation of the residuals.  See, for example, (5.4.4) on page 293 of
#  HM.  The function also assumes that the median of the raw residuals has
#  been used to estimate the intercept parameter.  In cases where the variances
#  may be non-positive, a least-squares type standardization is used.  However,
#  this is rarely needed.
#
#References:
#  Hettmansperger, T. P. and McKean, J. W. Robust nonparametric statistical methods.
#  Edward Arnold, London, 1998.
#
#  Naranjo, Joshua; McKean, Joseph W.; Sheather, Simon J.; Hettmansperger, Thomas P.
#  (1994). The use and interpretation of rank-based residuals.  J. Nonparametr. Statist.
#  3 no. 3-4, 323-341.
# 
#See Also:
#  mad, pairup, varcov.gr, varcov.hbr, wilcoxontau, taustar, studres.hbr
#
#Examples:
#> w1=wwfit(stars[,1],stars[,2],wilwts(stars[,1]))
#> w2=wwfit(stars[,1],stars[,2],grwts(stars[,1]))
#
#> # Wilcoxon (studentized) residuals
#> studres.gr(stars[,1],w1$weights,w1$resid)
# [1]  1.698464e-01  1.103038e+00 -3.742886e-01  1.103038e+00  5.654946e-02
# [6]  5.922285e-01 -1.138355e+00  3.810108e-01  6.138102e-01 -6.856994e-16
#[11]  3.138284e-01  5.542378e-01  5.454215e-01 -1.932241e+00 -1.386365e+00
#[16] -7.976069e-01 -1.925509e+00 -1.415705e+00 -1.554842e+00  5.773273e-01
#[21] -1.201132e+00 -1.448109e+00 -1.044846e+00 -3.294404e-01 -1.470635e-01
#[26] -6.739872e-01 -7.689204e-01 -3.323723e-01 -1.238059e+00  8.343941e-01
#[31] -1.073607e+00  1.095349e-01  2.135913e-01  1.236075e+00 -1.307730e+00
#[36]  9.641157e-01  8.722538e-02  2.135913e-01  2.112535e-01  7.397060e-01
#[41] -7.647593e-01 -3.381621e-02  4.366518e-01  3.991469e-01  7.848176e-01
#[46] -1.575200e-01 -9.212266e-01
#> plot(studres.gr(stars[,1],w1$weights,w1$resid))
#> abline(h=c(-2,2))
#
#> # GR (studentized) residuals
#> studres.gr(stars[,1],w2$weights,w2$resid)
# [1]  0.6832533  0.7330700  0.7032111  0.7330700  1.0482087  0.6608825
# [7]  2.6338108 -0.5586992  2.2454558  0.4257030  5.8043738  0.8230207
#[13]  0.4393625  0.2309870 -1.0930360 -1.1656892 -1.5141064 -2.1085260
#[19] -0.9245429  6.0953364 -0.8078962 -1.1880826 -1.5428239 -1.0100813
#[25]  0.1277571 -0.9771218 -0.1425699 -0.1533085 -0.3443044  6.3930775
#[31] -1.2775709 -0.9210367  0.1515091  6.8227429 -0.5315005  0.0000000
#[37] -0.6882846  0.1515091 -0.4877249  1.1066727 -0.8091282 -0.2294280
#[43]  0.1157052  0.4372119  0.2835305 -0.4198965 -1.3542565
#> plot(studres.gr(stars[,1],w2$weights,w2$resid))
#> abline(h=c(-2,2))
#
#Function:
#  studres.hbr
#
#Description:
#  This function calculates the (internally) studentized residuals
#  for the HBR estimate as discussed on page 323 of Hettmansperger and
#  McKean (1998). See also, Section 5.1.3 of Chang et al. (1994).
#  These residuals can be used to flag potential outliers (e.g. if the
#  absolute value of the studentized residual exceeds 2).
#
#Usage:
#  studres.hbr(x, bmat, res, delta=0.80, center=T)
#
#Arguments:
#       x: a n x p design matrix.
#    bmat: a n x n matrix that contains the weights, bij.  This matrix is
#          returned as the "weights" list component in the function "wwfit".
#     res: a n x 1 vector of residuals.  This vector is returned as the list
#          component "residuals" in the function "wwfit".
#   delta: represents the quantile used in the estimate of wilcoxontau.  The
#          default is 0.80.
#  center: should the design matrix be centered first.  The default is TRUE.
#  
#Value:
#  An n x 1 vector of studentized residuals defined as res/sqrt(diag(v)) where
#  v represents a first-order approximation of the variance-covariance matrix
#  of the raw (HBR) residuals.
#
#Details:
#  The variance-covariance residual matrix estimate is based on a first-order
#  approximation of the residuals.  See, for example, Section 5.1.3 of Chang 
#  et al. The function also assumes that the median of the raw residuals has
#  been used to estimate the intercept parameter.  In cases where the variances
#  may be non-positive, a least-squares type standardization is used.  However,
#  this is rarely needed.
#
#References:
#  Chang, W.H., McKean, J.W., Naranjo, J.D., and Sheather, S.J. (1999). "High
#  Breakdown Rank Regression". Journal of the American Statistical Association,
#  94(445), 205-219.
#
#  Hettmansperger, T. P. and McKean, J. W. Robust nonparametric statistical methods.
#  Edward Arnold, London, 1998.
# 
#See Also:
#  mad, pairup, varcov.gr, varcov.hbr, wilcoxontau, taustar, studres.gr
#
#Examples:
#> w1=wwfit(stars[,1],stars[,2],wilwts(stars[,1]))
#> w2=wwfit(stars[,1],stars[,2],grwts(stars[,1]))
#> w3=wwfit(stars[,1],stars[,2],hbrwts(stars[,1],stars[,2]))
#
#> # Wilcoxon (studentized) residuals
#> studres.gr(stars[,1],w1$weights,w1$resid)
# [1]  1.698464e-01  1.103038e+00 -3.742886e-01  1.103038e+00  5.654946e-02
# [6]  5.922285e-01 -1.138355e+00  3.810108e-01  6.138102e-01 -6.856994e-16
#[11]  3.138284e-01  5.542378e-01  5.454215e-01 -1.932241e+00 -1.386365e+00
#[16] -7.976069e-01 -1.925509e+00 -1.415705e+00 -1.554842e+00  5.773273e-01
#[21] -1.201132e+00 -1.448109e+00 -1.044846e+00 -3.294404e-01 -1.470635e-01
#[26] -6.739872e-01 -7.689204e-01 -3.323723e-01 -1.238059e+00  8.343941e-01
#[31] -1.073607e+00  1.095349e-01  2.135913e-01  1.236075e+00 -1.307730e+00
#[36]  9.641157e-01  8.722538e-02  2.135913e-01  2.112535e-01  7.397060e-01
#[41] -7.647593e-01 -3.381621e-02  4.366518e-01  3.991469e-01  7.848176e-01
#[46] -1.575200e-01 -9.212266e-01
#> plot(studres.gr(stars[,1],w1$weights,w1$resid))
#> abline(h=c(-2,2))
#
#> # GR (studentized) residuals
#> studres.gr(stars[,1],w2$weights,w2$resid)
# [1]  0.6832533  0.7330700  0.7032111  0.7330700  1.0482087  0.6608825
# [7]  2.6338108 -0.5586992  2.2454558  0.4257030  5.8043738  0.8230207
#[13]  0.4393625  0.2309870 -1.0930360 -1.1656892 -1.5141064 -2.1085260
#[19] -0.9245429  6.0953364 -0.8078962 -1.1880826 -1.5428239 -1.0100813
#[25]  0.1277571 -0.9771218 -0.1425699 -0.1533085 -0.3443044  6.3930775
#[31] -1.2775709 -0.9210367  0.1515091  6.8227429 -0.5315005  0.0000000
#[37] -0.6882846  0.1515091 -0.4877249  1.1066727 -0.8091282 -0.2294280
#[43]  0.1157052  0.4372119  0.2835305 -0.4198965 -1.3542565
#> plot(studres.gr(stars[,1],w2$weights,w2$resid))
#> abline(h=c(-2,2))
#
#> # HBR (studentized) residuals
#> studres.hbr(stars[,1],w3$weights,w3$resid)
# [1]  0.68141699  1.14819262  0.50340675  1.14819262  0.88453270  0.82974507
# [7]  3.72337389 -0.04838975  1.88044789  0.44983878  9.07013184  0.91702891
#[13]  0.67298687 -0.54215979 -1.04762774 -0.90040642 -1.52488679 -1.75062198
#[19] -0.99186343  9.35713657 -0.79224967 -1.13271260 -1.24055683 -0.63506446
#[25]  0.19927934 -0.73031510 -0.19816762 -0.05515384 -0.50336022  9.66221584
#[31] -1.06721936 -0.40664529  0.35192058 10.10858329 -0.64973229  0.62369290
#[37] -0.26075255  0.35192058 -0.07683403  1.17288945 -0.64607858  0.00000000
#[43]  0.42459843  0.60831384  0.70308785 -0.17363679 -1.07048540
#> plot(studres.hbr(stars[,1],w3$weights,w3$resid))
#> abline(h=c(-2,2))
#
#Function:
#   cellmntest
#
#Description:
#   Returns the (Wilcoxon) drop in dispersion test results and a robust
#   ANOVA table for cell means model.
#
#Usage:
#   cellmntest = function(y,levels,
#                         amat=cbind(rep(1,max(levels)-1),-1*diag(max(levels)-1)),
#                         delta=.80,param=2,print.tbl=T)
#
#Argurements:
#  y:         n x 1 vector of responses.
#  levels:    cell addresses (indices) for the responses.
#  amat:      q x p hypothesis matrix.  p is the number of cells.
#  delta:     bandwidth parameter for estimate of wilcoxontau.  Values should be
#             between .8 and .99.  Larger values result in a more conservative test.
#             If n/p is less than 5, recommend delta = .95.  Default is .8.
#  param:     a parameter for Huber's degrees of freedom correcton
#             as discussed in Huber (1981, p.174).  It is used in the estimator
#             of wilcoxontau.  Smaller values lead to more conservative estimates
#             of wilcoxontau.  The default is 2.  This is the benchmark for
#             labeling a standardized residual as a potential outlier.
#  print.tbl: if true, a robust ANOVA is printed, else results are returned.
#
#Value:
#
#  ANOVA     table of the hypothesis amat%*%mu = 0.
#  dred:     Value of dispersion function at reduced model (full model
#            constrained by amat%*%beta = 0).
#  dfull     Value of dispersion function at full model.
#  tauhat:   Estimate of scale parameter tau at the full model.
#  fr:       Wilcoxon drop in dispersion test.  To be compared with
#            F-critical values wiyth q and n-(p+1) degrees of freedom.
#  pval:     p-value using test statistic fr.
#
#Details:
#  Let the mu be the vector of cell means for the responses.   Let W
#  denote the incidence matrix corresponding to the vector of levels.
#  Then the full model is   Y = W%*%mu + e.   The hypotheses are
#           H_0: amat%*%mu = 0    versus  H_A: amat%*%mu not= 0.
#  The test statistic is
#                fr = [(dred-dfull)/q]/[tauhat/2]
#  The p-value is
#                pval = P[ F(q,n-(p+1)) >= fr],
#  where F(q,n-(p+1)) denotes a random variable with an F-distribution
#  with q and n-(p+1) degrees of freedom.   q*fr has an asymptotic chi^2
#  distribution under H_0 and a noncentral chi^2 under local alternatives.
#  See Section 3.6 and Chapter 4 of Hettmansperger and McKean (1998).
#
#References:
#  Hettmansperger, T. P. and McKean, J. W. Robust nonparametric statistical methods.
#    Edward Arnold, London, 1998.
#
#  Hollander, M. and Wolfe, D. A. Nonparametric statistical methods.
#    John Wiley & Sons, New York, 1999.
#
#See Also:
#  droptest, regrtest.
#
#Example:
#  This routine is easiest to use for a one-way design, but it can be used
#  for more complicated designs, also.  The example is based on a $2 \times 3$ 
#  crossed design.
#
#Consider the following data from a $2 \times 3$ crossed design:
#
#            B
#----------------------
#   (1)    (2)    (3) 
#A 53.7   51.3   50.2  
#  30.0   66.9   
#-----------------------
#   (4)    (5)    (6) 
#A  56.3   56.4   38.5
#   70.2   75.7 
#----------------------
#
#where the number in parantheses is the cell address.
#Suppose we are interested in testing for interaction.
#For a $2 \times 3$ this test has two degrees of freedom.
#We chose the two linearly independent
#contrasts given by:
#
#   (\mu_1 + \mu_5) - (\mu_2 + \mu_4)  =   0 
#   (\mu_2 + \mu_6) - (\mu_3 + \mu_5)  =   0 
#
#Hence, the marix of contrasts is given by
#1 -1  0  -1  1  0 
#0  1  -1  0  -1  1 
#
#> y
# [1] 53.7 51.3 50.2 30.0 66.9 56.3 56.4 38.5 70.2 75.7
#> levels
# [1] 1 2 3 1 2 4 5 6 5 6
#> amat
#   [,1] [,2] [,3] [,4] [,5] [,6]
#r1    1   -1    0   -1    1    0
#r2    0    1   -1    0   -1    1
#
#see = cellmntest(y,levels,amat)
#
#          RD DF     MRD     TS   PVAL
#H0    0.8033  2  0.4017 0.0124 0.9877
#Error     NA  5 32.4526     NA     NA
#
#Function:
#  pwcomp
#
#Description:   
#  This function conducts pairwise Wilcoxon-based drop in dispersion
#  tests for the one-way analysis of variance model.  It is typically
#  used ater cellmntest, if cellmntest yields a significant result.
#
#Usage:
#  pwcomp = function(y,levels,delta=.80,param=2)
#
#Arguements:
#  y:         n x 1 vector of responses.
#  levels:    cell addresses (indices) for the responses.
#  delta:     bandwidth parameter for estimate of wilcoxontau.  Values should be
#             between .8 and .99.  Larger values result in a more conservative test.
#             If n/p is less than 5, recommend delta = .95.  Default is .8.
#  param:     a parameter for Huber's degrees of freedom correcton
#             as discussed in Huber (1981, p.174).  It is used in the estimator
#             of wilcoxontau.  Smaller values lead to more conservative estimates
#             of wilcoxontau.  The default is 2.  This is the benchmark for
#             labeling a standardized residual as a potential outlier.
#
#Value:
#  A column vector which represents the p-value for the drop in dispersion test of
#  mu_i versus mu_j.  The order of the rows is given as (1,2), (1,3), ..., (1,p),
# (2,3), (2,4), ..., (2,p), ..., (p-1,p) where p represents the number of groups.
#
#Details:
#  The function simply calls cellmntest "p choose 2" times where each call uses
#  a corresponding amat argument.  Here, amat=(0,0,1,0,...,0,-1,0,...,0) where
#  the 1 corresponds to mu_i and the -1 corresponds to mu_j.  Note that this is
#  essentially equivalent to the protected LSD method.  See, for example, page 111
#  of Kuehl. 
#
#References:
#  Kuehl, Robert 0. Design of Experiments: Statistical Principles of Research Design
#    and Analysis.  Duxbury Press, Pacific Grove, CA, 2000.
#
#See Also:
#  cellmntest.
#
#Example:
#> y
# [1]  52  67  54  69 116  79  68  47 120  73  36  34  47 125  30  31  30  59  33
#[20]  98  52  55  66  50  58 176  91  66  61  63  62  71  41 118  48  82  65  72
#[39]  49
#> x
# [1] 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4
#[39] 4
#> t(pwcomp(y,x))
#      G1-G2  G1-G3  G1-G4  G2-G3  G2-G4  G3-G4
#PVAL 0.0031 0.5984 0.5433 0.0131 0.0173 0.8472
#
#Function:
#   cellmnxy
#
#Description:
#   Obtains the incidence matrix for the function cellmntest.
#
#Usage:
#   cellmnxy (levels)
#
#Arguements:
#   levels:    Cell addresses for the corresponding vector of responses.
#
#Value:
#   cellmnxy:  The incidence matrix.
#
#See Also:
#   cellmntest
#
#Example:
#
#> y
# [1] 53.7 51.3 50.2 30.0 66.9 56.3 56.4 38.5 70.2 75.7
#> levels
# [1] 1 2 3 1 2 4 5 6 5 6
#>
#> cellmnxy(y,levels)
#      [,1] [,2] [,3] [,4] [,5] [,6]
# [1,]    1    0    0    0    0    0
# [2,]    0    1    0    0    0    0
# [3,]    0    0    1    0    0    0
# [4,]    1    0    0    0    0    0
# [5,]    0    1    0    0    0    0
# [6,]    0    0    0    1    0    0
# [7,]    0    0    0    0    1    0
# [8,]    0    0    0    0    0    1
# [9,]    0    0    0    0    1    0
#[10,]    0    0    0    0    0    1
#
#Function:
#  droptest
#
#Description:
#  The function droptest obtains the drop in dispersion test
#  for a general linear hypothesis.
#
#Usage:
#  droptest(xmat,y,amat,delta=.80,param=2,print.tbl=T)
#
#Arguements:
#  xmat:      n x p design matrix with no intercept column.  Assumed
#             to have full column rank.
#  y:         n x 1 vector of responses.
#  amat:      q x p hypothesis matrix of full row rank.  
#  delta:     bandwidth parameter for estimate of wilcoxontau.  Values should be
#             between .8 and .99.  Larger values result in a more conservative test.
#             If n/p is less than 5, recommend delta = .95.  Default is .8.
#  param:     a parameter for Huber's degrees of freedom correcton
#             as discussed in Huber (1981, p.174).  It is used in the estimator
#             of wilcoxontau.  Smaller values lead to more conservative estimates
#             of wilcoxontau.  The default is 2.  This is the benchmark for
#             labeling a standardized residual as a potential outlier.
#  print.tbl: if true, a robust ANOVA is printed.
#
#
#Value:
#  list(list(full=full,dred=dred,dfull=dfull,tauhat=tauhat,q=q,fr=fr,pval=pval))
#
#  dred:     Value of dispersion function at reduced model (full model
#            constrained by amat%*%beta = 0).
#  dfull     Value of dispersion function at full model.
#  tauhat:   Estimate of scale parameter tau at the full model.
#  q:        Number of constraints (no. of rows of amat).
#  fr:       Wilcoxon drop in dispersion test.  To be compared with
#            F-critical values wiyth q and n-(p+1) degrees of freedom.
#  pval:     p-value using test statistic fr.
#
#Details:
#  Let the full model be   Y = xmat%*%beta + e.   The hypotheses are
#           H_0: amat%*%beta = 0    versus  H_A: amat%*%beta not= 0.
#  The test statistic is
#                fr = [(dred-dfull)/q]/[tauhat/2]
#  The p-value is
#                pval = P[ F(q,n-(p+1)) >= fr],
#  where F(q,n-(p+1)) denotes a random variable with an F-distribution
#  with q and n-(p+1) degrees of freedom.   q*fr has an asymptotic chi^2
#  distribution under H_0 and a noncentral chi^2 under local alternatives.
#  See Section 3.6 of Hettmansperger and Mckean (1998).
#
#References:
#  Hettmansperger, T. P. and McKean, J. W. Robust nonparametric statistical methods.
#    Edward Arnold, London, 1998.
#
#  Hollander, M. and Wolfe, D. A. Nonparametric statistical methods.
#    John Wiley & Sons, New York, 1999.
#
#See Also:
#  cellmntest, regrtest, redmod
#
#Examples:
#
#Consider the following one-way ANOVA data set taken from page 236
#of Hettmansperger and McKean (1998).  Thirty-nine quail were randomly 
#assigned to four diets, each diet containing a different
#drug compound.  The drug compounds are labeled: I, II, III, and IV.
#At the end of the prescribed experimental time the LDL cholesterol of each quail
#was measured.
#
#Data:
#I   52  67  54  69  116  79  68  47  120  73
#II  36  34  47  125  30  31  30  59  33  98 
#III 52  55  66  50  58  176  91  66  61  63
#IV  62  71  41  118  48  82  65  72  49  
#
#As a design matrix miuns the intercept column was composed of
#three coulmns. The first column is the incidence vector for 
#Level II, the second for Level III, and the third for Level IV.
#Our hypothesis of interest is that there is no level effect; hence,
#the hypothesis matrix is the identity matrix.
#
#Based on the ouput shown below, the test statistic is 3.844.
#This should be compared with an $F$-distribution with 3 and 35
#degrees of freedom.
#
#> amat
#   [,1] [,2] [,3]
#r1    1    0    0
#r2    0    1    0
#r3    0    0    1
#
#> temp = droptest(xmat,y,amat)^C
#> temp=droptest(x,y,amat)
#
#
#           RD DF     MRD     TS   PVAL
#H0    108.611  3 36.2037 3.7896 0.0185
#Error      NA 36  9.5535     NA     NA
#
#Function:
#   fitdiag
#
#Description:
#   Preforms diagnostics between two robust fits.
#
#Usage:
#   fitdiag (x,y,est=c("WIL","GR"),delta=0.80,param=2,conf=0.95)
#
#Arguments:
#   x:      Design matrix (need not be centered).
#   y:      Vector of responses.
#   est:    vector containing any two of "LS", "WIL", "GR", or "HBR".
#   delta:  parameter for estimate of tau
#   param:  Huber's correction for df correction for estimate of tau
#   conf:   Confidence coefficient used in CI estimate of taustar
#
#Value:
#   tdbeta: total difference in fit (standardized at WIL)
#   bmtd:   Benchmark for tdbeta
#   cfit:   Caswise difference diagnostics
#   bmcf:   Benchmark for cfit
#
#Details:
#   This function returns the diganostics (TDBETA and CFIT) between
#   any two comparisons between LS, WIL, GR, and HBR.  All comparisons
#   are standardized using the WIL fit.
#   See, for example, McKean, Naranjo and Sheather in the articles: 
#   (1996, Computational Stat., 223-243) and
#   (1996,Comm. in Stat-Theory, 2575-2595) or section 5.5 of HM.
#
#References:
#   Hettmansperger, T. P. and McKean, J. W. Robust nonparametric statistical methods.
#     Edward Arnold, London, 1998.
#
#   McKean, J. W., Naranjo, J. D. and Sheather, S. J. (1996a), Diagnostics
#     to detect differences in robust fits of linear models, Computational
#     Statistics, 11, 223-243.
#
#   McKean, J. W., Naranjo, J. D. and Sheather, S. J. (1996b),
#     An efficient and high breakdown procedure for model criticism, 
#     Communications in Statistics, Part A-Theory and Methods, 25, 2575-2595.
#
#   McKean, J. W., Naranjo, J. D. and Sheather, S. J. (1999),
#      Diagnostics for comparing robust and least squares fits,
#      Journal of Nonparametric Statistics, 11, 161-188.
#
#See Also:
#   plotfitdiag
#
#Examples:
#   Example 1: Difference between HBR and Wilcoxon fits.
#
#> x
#      c1 c2
# [1,] 40 63
# [2,] 46 43
# [3,] 42 46
# [4,] 47 56
# [5,] 72 50
# [6,] 56 46
# [7,] 63 56
# [8,] 37 68
# [9,] 35 41
#[10,] 82 92
#> y
# [1] 897 709 765 902 922 794 941 944 671 813
#
#> fitdiag(x,y,est=c("WIL","HBR"))
#$tdbeta
#[1] 9.800308
#
#$bmtd
#[1] 3.6
#
#$cfit
# [1] -0.5765829  1.0326939  0.9544706 -0.3337067 -0.5998241  0.2763236
# [7] -1.0327780 -0.7327053  1.5251840 -3.0071011
#
#$bmcf
#[1] 1.095445
#
#   Example 2:  Difference between Wilcoxon and LS fits.
#
#> x
#      c1 c2
# [1,] 40 63
# [2,] 46 43
# [3,] 42 46
# [4,] 47 56
# [5,] 72 50
# [6,] 56 46
# [7,] 63 56
# [8,] 37 68
# [9,] 35 41
#[10,] 82 92
#> y
# [1] 897 709 765 902 922 794 941 944 671 813
#
#> fitdiag(x,y,est=c("WIL","LS"))
#$tdbeta
#[1] 0.07609967
#
#$bmtd
#[1] 3.6
#
#$cfit
# [1] 0.09498017 0.24395623 0.21596049 0.20834416 0.24355930 0.27049185
# [7] 0.26614446 0.04858123 0.17560031 0.09359674
#
#$bmcf
#[1] 1.095445
#
#Function:
#   plotfitdiag
#
#Description:
#   Obtains casewise plot of change in fits based on the results of the function plotfitdiag.
#
#Usage:
#   plotfitdiag (result)
#   where
#   result = fitdiag(x,y,est=c("WIL","GR"),delta=0.80,param=2,conf=0.95)
#
#Arguments:
#   result: Results of the function plotfitdiag.
#
#Value:
#   Casewise plot of change in fits, including display of total change in fits.
#
#See Also:
#   fitdiag
#
#Function:
#   redmod
#
#Description:
#   Obtains a reduced model design matrix used in the function droptest.
#   The function only needs a matrix whose column space is a basis matrix
#   for the reduced model space.
#
#Usage:
#  redmod (xmat,amat)
#
#Arguements:
#   xmat:    n x p design matrix.
#   amat:    q x p hypothesis matrix.
#
#Value:
#   redmod:   n x (p-q) reduced model design matrix.
#
#Details:
#   The function utilizes results from QR decompositions;
#   see Theorem 3.7.2 of Hettmansperger and McKean (1998).
#
#References:
#   Hettmansperger, T. P. and McKean, J. W. Robust nonparametric statistical methods.
#   Edward Arnold, London, 1998.
#
#Example:
#> xmat
#     [,1] [,2] [,3]
#[1,]  109   70   96
#[2,]  108  104  101
#[3,]   98   80   96
#[4,]   99   73  108
#[5,]  115   88   95
#[6,]  132   85   95
#> r1=c(1,0,0)
#> r2=c(0,0,1)
#> amat=rbind(r1,r2)
#> redmod(xmat,amat)
#     [,1]
#[1,]  -70
#[2,] -104
#[3,]  -80
#[4,]  -73
#[5,]  -88
#[6,]  -85
#
#Function:
#   regrtest
#
#Description:
#   Wilcoxon drop in dispersion test of regression significance.
#
#Usage:
#   regrtest (xmat,y,delta=0.80,param=2,print.tbl=T)
#
#Arguements:
#   xmat:     n x p design matrix with no intercept column.  Assumed
#             to have full column rank.
#   y:        n x 1 vector of responses.
#   delta:    bandwidth parameter for estimate of wilcoxontau.  Values should be
#             between .8 and .99.  Larger values result in a more conservative test.
#             If n/p is less than 5, recommend delta = .95.  Default is .8.
#  param:     a parameter for Huber's degrees of freedom correcton
#             as discussed in Huber (1981, p.174).  It is used in the estimator
#             of wilcoxontau.  Smaller values lead to more conservative estimates
#             of wilcoxontau.  The default is 2.  This is the benchmark for
#             labeling a standardized residual as a potential outlier.
#  print.tbl  if true, a robust ANOVA is printed.  If False test inoformation is
#             returned.
#
#Value:
#   ANOVA     Table is printed out for Wilcoxon droptest of regression significance.
#   dred:     Value of dispersion function at reduced model (full model
#             constrained by amat%*%beta = 0).
#   dfull     Value of dispersion function at full model.
#   tauhat:   Estimate of scale parameter tau at the full model.
#   fr:       Wilcoxon drop in dispersion test.  To be compared with
#             F-critical values wiyth q and n-(p+1) degrees of freedom.
#   pval:     p-value using test statistic fr.
#
#Details:
#  Let the full model be   Y = xmat%*%beta + e.   The hypotheses are
#           H_0: beta = 0    versus  H_A: beta not= 0.
#  The test statistic is
#                fr = [(dred-dfull)/p]/[tauhat/2]
#  The p-value is
#                pval = P[ F(p,n-(p+1)) >= fr],
#  where F(q,n-(p+1)) denotes a random variable with an F-distribution
#  with q and n-(p+1) degrees of freedom.   p*fr has an asymptotic
#  chi^2 distribution under H_0 and a noncentral chi^2 under local
#  alternatives.  See Section 3.6 of Hettmansperger and McKean (1998).
#
#References:
#  Hettmansperger, T. P. and McKean, J. W. Robust nonparametric
#    statistical methods.  Edward Arnold, London, 1998.
#
#  Hollander, M. and Wolfe, D. A. Nonparametric statistical methods.
#    John Wiley & Sons, New York, 1999.
#
#See Also:
#  cellmntest, droptest, redmod
#
#Example:
#
#> x=matrix(rnorm(20),ncol=2)
#> y = rnorm(10)
#> regrtest(x,y)
#
#         RD DF    MRD    TS   PVAL
#H0    0.649  2 0.3245 1.002 0.4143
#Error        7 0.3238             
#
#Function:
#   wilcoxonpseudo
#
#Description:
#   Obtains the Wilcoxon pseudo-observations.
#
#Usage:
#  wilcoxonpseudo (x,y,delta=.80,param=2)
#
#Arguements:
#   x:    n x p design matrix.
#   y:    n x 1 vector of responses.
#????????????????????????
#Value:
#   redmod   n x (p-q) reduced model design matrix.
#
#Details:
#   The function utilizes results from QR decompositions;
#   see Theorem 3.7.2 of Hettmansperger and McKean (1998).
#
#References:
#   Hettmansperger, T. P. and McKean, J. W. Robust nonparametric statistical methods.
#   Edward Arnold, London, 1998.
#
#Example:
#> xmat
#     [,1] [,2] [,3]
#[1,]  109   70   96
#[2,]  108  104  101
#[3,]   98   80   96
#[4,]   99   73  108
#[5,]  115   88   95
#[6,]  132   85   95
#> r1=c(1,0,0)
#> r2=c(0,0,1)
#> amat=rbind(r1,r2)
#> redmod(xmat,amat)
#     [,1]
#[1,]  -70
#[2,] -104
#[3,]  -80
#[4,]  -73
#[5,]  -88
#[6,]  -85
#
#Function:
#   wildisp
#
#Description:
#   Dispersion function based on Wilcoxon scores.
#
#Usage:
#   wildisp(resid)
#
#Argument:
#   resid:     Vector at which the dispersion function is evaluated.
#
#Value:
#   wildisp:   Value of dispersion function at resid.
#
#Details:
#   The function uses standardized scores, so it corresponds
#   to the dispersion function in formula (3.2.6) given on page 146
#   of Hettmansperger and McKen (1998) for Wilcoxon scores.
#
#References:
#   Hettmansperger, T. P. and McKean, J. W. Robust nonparametric statistical methods.
#     Edward Arnold, London, 1998.
#
#Example:
#
#> y = rnorm(10)
#> y
# [1]  1.2540313 -0.5230302  0.7303705 -1.8067392  1.3255736 -1.2026831
# [7] -1.7971791 -0.7946445 -0.4380870 -0.4843025
#>
#> wildisp(y)
#[1] 11.01341
