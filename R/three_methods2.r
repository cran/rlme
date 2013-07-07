GEER_est2 <-
function (x, y, I, sec, mat, school, section = 1, weight = "wil", 
    rprpair = "hl-disp") 
{
    weight = tolower(weight)
    if (weight == "wil") {
        weight = 1
    }
    if (weight == "hbr") {
        weight = 2
    }
    fitw = wwest(x, y, print.tbl = F)
    b0 = fitw$tmp1$coef
    ehat0 = y - cbind(1, x) %*% b0
    fitvc = rprmeddis2(I, sec, mat, ehat = ehat0, location, scale, 
        rprpair = rprpair)
    sigmaa2 = fitvc$siga2
    sigmae2 = fitvc$sigmae2
    x = cbind(1, x)
    collb <- b0
    collsigma <- c(sigmaa2, sigmae2)
    iter <- 0
    chk <- 1
    b2 <- b0
    max.iter <- 2
    ww <- 0
    while ((chk > 1e-04) && (iter < max.iter)) {
        iter <- iter + 1
        b1 <- b2
        n = sum(mat)
        sigmay = sigymake2(I, sec, mat, sigmaa2, sigmae2)
        sigma12inv = matrix(sigmay$sigy12i, ncol = n)
        siggma12 <- matrix(sigma12(sigmay$sigy2), ncol = n)
        ystar = sigma12inv %*% y
        xstar = sigma12inv %*% x
        ehats = ystar - xstar %*% b2
        med <- median(ehats)
        ahats = scorewil(ehats)$scorewil
        yss = y - siggma12 %*% rep(1, n) * med
        ys = ystar - siggma12 %*% rep(1, n) * med
        if (weight == 1) {
            w <- weightf(ehats, ahats, med)$w
        }
        if (weight == 2) {
            w <- hbrwts_gr(xstar, ehats)
        }
        ww <- rbind(ww, t(w))
        w <- diag(as.vector(w))
        r <- sigma12inv %*% w %*% sigma12inv
        b2 <- solve(t(x) %*% r %*% x, tol = tol) %*% (t(x) %*% 
            r %*% yss)
        ehat = y - x %*% b2
        fitvc <- rprmeddis2(I, sec, mat, ehat = ehat, location, 
            scale, rprpair = rprpair)
        sigmaa2 = fitvc$siga2
        sigmae2 = fitvc$sigmae2
        chk <- sum((b2 - b1)^2)/sum(b1^2)
        collb <- rbind(collb, c(b2))
        collsigma <- rbind(collsigma, c(sigmaa2, sigmae2))
    }
    iter = iter + 1
    b <- collb[iter, ]
    theta = b
    ehat = y - x %*% b
    ahat = scorewil(ehat)$scorewil
    tauhat = wilcoxontau(ehat, p = dim(x)[2])
    fitvc <- rprmeddis2(I, sec, mat, ehat = ehat, location, scale, 
        rprpair = rprpair)
    sigmaa2 = fitvc$siga2
    sigmae2 = fitvc$sigmae2
    sigma <- c(sigmaa2, sigmae2)
    effect_sch = fitvc$frei
    effect_err = fitvc$free
    sigmay = sigymake2(I, sec, mat, sigmaa2, sigmae2)
    sigma12inv = matrix(sigmay$sigy12i, ncol = n)
    siggma12 <- matrix(sigma12(sigmay$sigy2), ncol = n)
    if (weight == 1) {
        w <- weightf(ehats, ahats, med)$w
    }
    if (weight == 2) {
        w <- hbrwts_gr(xstar, ehats)
    }
    w <- diag(as.vector(w))
    r <- sigma12inv %*% w %*% sigma12inv
    ystar = sigma12inv %*% y
    xstar = sigma12inv %*% x
    ehats = ystar - xstar %*% b
    ahats = scorewil(ehats)$scorewil
    tauhats = wilcoxontau(ehats, p = dim(x)[2])
    m01 = solve(t(x) %*% r %*% x, tol = tol)
    m03 = solve(t(x) %*% sigma12inv %*% sigma12inv %*% x, tol = tol)
    sisi3 = diag(rep(1, n))
    m15 = (t(x) %*% sigma12inv %*% sisi3 %*% sigma12inv %*% x)
    varb = tauhats^2 * m03 %*% m15 %*% m03
    se31 = sqrt(diag(varb))
    pp <- dim(x)[2] + 1
    cc = jrfit2(pp, ehats, ahats, block = school)
    rho1 = cc$rhos
    v1 = 1
    v2 = rho1
    sisi4 = 0
    ublock <- unique(school)
    for (i in ublock) {
        sisi4 = adiag(sisi4, Bmat2(v1, v2, school[school == i]))
    }
    sisi4 = sisi4[2:dim(sisi4)[1], 2:dim(sisi4)[2]]
    sisi5 = matrix(sisi4, ncol = n)
    m16 = (t(x) %*% sigma12inv %*% sisi4 %*% sigma12inv %*% x)
    if (weight == 2) {
        varb = m01 %*% m16 %*% m01
        se41 = sqrt(diag(varb))
    }
    if (weight == 1) {
        varb = tauhats^2 * m03 %*% m16 %*% m03
        se41 = sqrt(diag(varb))
    }
    list(theta = theta, ses_AP = se31, ses_CS = se41, varb = varb, 
        sigma = sigma, ehat = ehat, effect_sch = effect_sch, 
        effect_err = effect_err, iter = iter, w = diag(w))
}

GR_est2 <-
function (x, y, I, sec, mat, school, section = 1, rprpair = "hl-disp") 
{
    init = T
    sigmaa2 = 1
    sigmae2 = 1
    thetaold = c(0)
    numstp = 50
    eps = 1e-04
    iflag2 = 0
    i = 0
    is = 0
    J = sum(sec)
    n = length(y)
    if (is.null(dim(x)[2])) {
        nopar_scale = dim(x)[2] + 1 + 2
    }
    else {
        nopar_scale = dim(x)[2] + 1 + 2
    }
    coll = matrix(rep(0, nopar_scale), ncol = nopar_scale)
    collresch = matrix(rep(0, I), ncol = I)
    collresec = matrix(rep(0, J), ncol = J)
    collreserr = matrix(rep(0, n), ncol = n)
    while (is == 0) {
        i = i + 1
        if (i == 1) {
            wilfit = wilstep2(I, sec, mat, init, y, x, sigmaa2, 
                sigmae2, thetaold, eps, iflag2, rprpair = rprpair)
            coll = rbind(coll, c(wilfit$theta, wilfit$sigmaa2, 
                wilfit$sigmae2))
            effect_sch = wilfit$rea
            effect_err = wilfit$ree
            thetaold = wilfit$theta
            ehat = wilfit$ehat
            theta = coll[dim(coll)[1], 1:(dim(x)[2] + 1)]
            sigma = c(wilfit$sigmaa2, wilfit$sigmae2)
        }
        sigmaa2 = wilfit$sigmaa2
        sigmae2 = wilfit$sigmae2
        init = F
        wilfit = wilstep2(I, sec, mat, init, y, x, sigmaa2, sigmae2, 
            thetaold, eps, iflag2, rprpair = rprpair)
        if ((wilfit$iflag2 == 1) | (i > numstp)) {
            is = 1
        }
        coll = rbind(coll, c(wilfit$theta, wilfit$sigmaa2, wilfit$sigmae2))
        effect_sch = wilfit$rea
        effect_err = wilfit$ree
        thetaold = wilfit$theta
        ehat = wilfit$ehat
        theta = coll[dim(coll)[1], 1:(dim(x)[2] + 1)]
        sigma = c(wilfit$sigmaa2, wilfit$sigmae2)
        i
    }
    sigmay = sigymake2(I, sec, mat, sigma[1], sigma[2])
    sigma12inv = matrix(sigmay$sigy12i, ncol = length(y))
    xstar = sigma12inv %*% cbind(1, x)
    ystar = sigma12inv %*% y
    ehats = ystar - xstar %*% theta
    ahats = scorewil(ehats)$scorewil
    ehat = y - cbind(1, x) %*% theta
    taus = taustar(ehats, p = dim(x)[2], conf = 0.95)
    tauhat = wilcoxontau(ehats, p = dim(x)[2])
    XXinv <- solve(crossprod(xstar))
    x_c = xstar - apply(xstar, 2, mean)
    aa = taus^2 * (XXinv %*% t(xstar)) %*% {
        rep(1, n) %*% solve(t(rep(1, n)) %*% rep(1, n)) %*% rep(1, 
            n)
    } %*% (xstar %*% XXinv)
    bb = tauhat^2 * XXinv %*% t(xstar) %*% {
        x_c %*% solve(t(x_c) %*% x_c) %*% t(x_c)
    } %*% (xstar %*% XXinv)
    varb = aa + bb
    ses <- sqrt(diag(varb))
    list(theta = theta, ses = ses, sigma = sigma, varb = varb, 
        ehat = ehat, ehats = ehats, effect_sch = effect_sch, 
        effect_err = effect_err, iter = i, coll = coll, xstar = xstar, 
        ystar = ystar)
}

JR_est2 <-
function (x, y, I, sec, mat, school, section = 1, rprpair = "hl-disp") 
{
    pp <- dim(x)[2] + 1
    wilfit = wilonestep(y, x)
    theta = wilfit$theta[1:(pp)]
    ehat = wilfit$ehat
    ahat = scorewil(ehat)$scorewil
    taus = wilfit$taus
    tauhat = wilfit$tauhat
    scale.fit = rprmeddis2(I, sec, mat, ehat, location, scale, 
        rprpair = rprpair)
    sigma = c(scale.fit$siga2, scale.fit$sigmae2)
    effect_sch = scale.fit$frei
    effect_err = scale.fit$free
    pp <- dim(x)[2] + 1
    cc = jrfit2(pp, ehat, ahat, block = school)
    var_alpha = taus^2 * (1/length(ehat)) * cc$sigmastar
    rho1 = cc$rhos
    v1 = 1
    v2 = rho1
    V = beta_var2(x, school, tauhat, v1, v2)$var
    ses <- c(sqrt(var_alpha), sqrt(diag(V)))
    theta <- as.vector(theta)
    sigma <- sigma
    list(theta = theta, ses = ses, varb = V, sigma = sigma, ehat = ehat, 
        effect_sch = effect_sch, effect_err = effect_err)
}
LM_est2 <-
function (x, y, dat, method = "REML") 
{
    model = as.formula(paste("y ~ 1 + ", paste(colnames(x), collapse = " + ")))
    fit.lme = lme(model, data = dat, random = ~1 | school, method = method)
    summary(fit.lme)
    theta0 <- extract.lme.cov2(fit.lme, dat, start.level = 1)$V[[1]][1, 
        2]
    theta1 <- extract.lme.cov2(fit.lme, dat, start.level = 1)$V[[1]][1, 
        1] - (theta0)
    sigma.l <- c(theta0, theta1)
    theta <- as.vector(summary(fit.lme)$tTable[, 1])
    intra_sch.lm <- (theta1)/(sum(sigma.l))
    ses <- as.vector(summary(fit.lme)$tTable[, 2])
    ehat <- as.vector(fit.lme$residuals[, 1])
    varb = fit.lme$varFix
    effect_sch = random.effects(fit.lme, level = 1)[, 1]
    effect_err = as.vector(fit.lme$residuals[, 2])
    list(theta = theta, ses = ses, varb = varb, sigma = sigma.l, 
        ehat = ehat, effect_sch = effect_sch, effect_err = effect_err)
}
