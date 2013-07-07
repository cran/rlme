stanresidgr <-
function (x, y, resid, delta = 0.8, param = 2, conf = 0.95) 
{
    xc = as.matrix(centerx(x))
    n = length(y)
    p = length(xc[1, ])
    pp1 = p + 1
    hc = diag(xc %*% ginv(t(xc) %*% xc) %*% t(xc))
    tau = wilcoxontau(resid, p, delta = 0.8, param = 2)
    taus = taustar(resid, p, conf = 0.95)
    deltas = sum(abs(resid))/(n - pp1)
    delta = wildisp(resid)/(n - pp1)
    sig = mad(resid)
    k1 = (taus^2/sig^2) * (((2 * deltas)/taus) - 1)
    k2 = (tau^2/sig^2) * (((2 * delta)/tau) - 1)
    s1 = sig^2 * (1 - (k1/n) - k2 * hc)
    s2 = s1
    s2[s1 <= 0] = sig^2 * (1 - (1/n) - hc[s1 <= 0])
    ind = rep(0, n)
    ind[s1 <= 0] = 1
    stanresid = resid/sqrt(s2)
    list(stanr = stanresid, ind = ind, rawresids = resid, tau = tau, 
        taustar = taus)
}
