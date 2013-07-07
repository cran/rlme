hbrwts_gr <-
  function (xmat, y, percent = 0.95,
            intest = ltsreg(xmat,y)$coef) 
  {
    # modified hbrwts function in ww and rlme.v2
    # need to install two packs:
    # install.packages("robustbase")
    # install.packages("MASS")  
    # library(robustbase)
    # library(MASS)
    # xmat= x, y is dependent
    # check it with older ww function
    
    xmat = as.matrix(xmat)
    
    aa=covMcd(xmat)

    if(dim(xmat)[2]==1){  
      robdis2=mahalanobis(xmat, aa$raw.center, aa$raw.cov)
    }  else
    {robdis2=aa$raw.mah}
    
    
    y = as.matrix(y)
    n = dim(xmat)[1]
    p = dim(xmat)[2]
    cut = qchisq(percent, p)
    resids = y - intest[1] - xmat %*% as.matrix(intest[2:(p + 1)])
    sigma = mad(resids)
    m = psi(cut/robdis2)
    a = resids/(sigma * m)
    c = (median(a) + 3 * mad(a))^2
    h = sqrt(c)/a
    ans = psi(abs(h))
    #    ans=psi(abs(h)^sqrt(sum(h^2)))
    ans
  }