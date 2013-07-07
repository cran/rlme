wilstep2 <-
function (I, sec, mat, init = F, y, x, sigmaa2 = 1, sigmae2 = 1, 
    thetaold = c(0), eps = 1e-04, iflag2 = 0, rprpair = "hl-disp") 
{
    location = scale = 2
    if (rprpair == "med-mad") {
        location = scale = 1
    }
    n = length(y)
    if (init == T) {
        fitw = wwest(x, y, print.tbl = F)
        theta = fitw$tmp1$coef
        ehat = fitw$tmp1$residuals
        fitvc = rprmeddis2(I, sec, mat, ehat, location, scale, 
            rprpair = rprpair)
        sigmaa2 = fitvc$siga2
        sigmae2 = fitvc$sigmae2
        rea = fitvc$frei
        ree = fitvc$free
    }
    else {
        sigmay = sigymake2(I, sec, mat, sigmaa2, sigmae2)
        sigma12inv = matrix(sigmay$sigy12i, ncol = n)
        ystar = sigma12inv %*% y
        xstar = sigma12inv %*% cbind(rep(1, n), x)
        fitw = wwest(xstar, ystar, print.tbl = F)
        yhat = ystar - fitw$tmp1$residuals
        fitcsp = projcsp(xstar, yhat)
        theta = fitcsp$betahat
        ehat = y - cbind(rep(1, n), x) %*% theta
        fitvc = rprmeddis2(I, sec, mat, ehat, location, scale, 
            rprpair = rprpair)
        sigmaa2 = fitvc$siga2
        sigmae2 = fitvc$sigmae2
        rea = fitvc$frei
        ree = fitvc$free
        chk = sum((theta - thetaold)^2)/sum(thetaold^2)
        if (chk < eps) {
            iflag2 = 1
        }
    }
    list(theta = theta, ehat = ehat, sigmaa2 = sigmaa2, sigmae2 = sigmae2, 
        rea = rea, ree = ree, iflag2 = iflag2)
}
