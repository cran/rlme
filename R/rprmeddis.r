#rpr <-
#function (ehat, school, section, rprpair = "hl-disp") 
#{
#    I = length(unique(factor(section)))
#    sec = as.vector(sec_vec(school, section))
#    mat = mat_vec(school, section)
#    rprpair = tolower(rprpair)
#    location = scale = 2
#    if (rprpair == "med-mad") {
#        location = scale = 1
#    }
#    return(rprmeddis(I, sec, mat, ehat, location, scale))
#}

rprmeddis <-
function (I, sec, mat, ehat, location, scale, rprpair = "hl-disp") 
{
    rprpair = tolower(rprpair)
    if (rprpair == "hl-disp") {
        location = scale = 2
    }
    if (rprpair == "med-mad") {
        location = scale = 1
    }
    if (location == 1) {
        matre = matrix(c(0, 0), ncol = 2)
        uhatij = 0 * mat
        ni = apply(mat, 1, sum)
        reij = 0 * mat
        rei = rep(0, I)
        rew = rep(0, sum(sec))
        frei = rep(0, I)
        uhati = rep(0, I)
        ehat2 = ehat
        ehat3 = 0 * ehat
        ic = 0
        id = 0
        ia = 0
        for (i in 1:I) {
            for (j in 1:sec[i]) {
                ehattmp = ehat[(ic + 1):(ic + mat[i, j])]
                ic = ic + mat[i, j]
                uhatij[i, j] = median(ehattmp)
            }
            uhati[i] = median(uhatij[i, 1:sec[i]])
            for (j in 1:sec[i]) {
                reij[i, j] = uhatij[i, j] - uhati[i]
            }
            for (j in 1:sec[i]) {
                ehat2[(id + 1):(id + mat[i, j])] = ehat[(id + 
                  1):(id + mat[i, j])] - reij[i, j]
                id = id + mat[i, j]
            }
            rei[i] = median(ehat2[(ia + 1):(ia + ni[i])])
            ia = ia + ni[i]
        }
        frei = rei - median(rei)
        ib = 0
        for (i in 1:I) {
            for (j in 1:sec[i]) {
                ib = ib + 1
                rew[ib] = reij[i, j]
            }
        }
        frew = rew - median(rew)
    }
    if (location == 2) {
        matre = matrix(c(0, 0), ncol = 2)
        uhatij = 0 * mat
        ni = apply(mat, 1, sum)
        reij = 0 * mat
        rei = rep(0, I)
        rew = rep(0, sum(sec))
        frei = rep(0, I)
        uhati = rep(0, I)
        ehat2 = ehat
        ehat3 = 0 * ehat
        ic = 0
        id = 0
        ia = 0
        for (i in 1:I) {
            for (j in 1:sec[i]) {
                ehattmp = ehat[(ic + 1):(ic + mat[i, j])]
                ic = ic + mat[i, j]
                if (length(ehattmp) == 1) {
                  uhatij[i, j] = ehattmp
                }
                else {
                  uhatij[i, j] = onesampwil(ehattmp, maktable = F, 
                    plotb = F)$est
                }
            }
            if (length(uhatij[i, 1:sec[i]]) == 1) {
                uhati[i] = uhatij[i, 1:sec[i]]
            }
            else {
                uhati[i] = onesampwil(uhatij[i, 1:sec[i]], maktable = F, 
                  plotb = F)$est
            }
            for (j in 1:sec[i]) {
                reij[i, j] = uhatij[i, j] - uhati[i]
            }
            for (j in 1:sec[i]) {
                ehat2[(id + 1):(id + mat[i, j])] = ehat[(id + 
                  1):(id + mat[i, j])] - reij[i, j]
                id = id + mat[i, j]
            }
            if (length(ehat2[(ia + 1):(ia + ni[i])]) == 1) {
                rei[i] = ehat2[(ia + 1):(ia + ni[i])]
            }
            else {
                rei[i] = onesampwil(ehat2[(ia + 1):(ia + ni[i])], 
                  maktable = F, plotb = F)$est
            }
            ia = ia + ni[i]
        }
        frei = rei - onesampwil(rei, maktable = F, plotb = F)$est
        ib = 0
        for (i in 1:I) {
            for (j in 1:sec[i]) {
                ib = ib + 1
                rew[ib] = reij[i, j]
            }
        }
        if (length(rew) == 1) {
            frew = rew - rew
        }
        else {
            frew = rew - onesampwil(rew, maktable = F, plotb = F)$est
        }
    }
    if (location == 3) {
        matre = matrix(c(0, 0), ncol = 2)
        uhatij = 0 * mat
        ni = apply(mat, 1, sum)
        reij = 0 * mat
        rei = rep(0, I)
        rew = rep(0, sum(sec))
        frei = rep(0, I)
        uhati = rep(0, I)
        ehat2 = ehat
        ehat3 = 0 * ehat
        ic = 0
        id = 0
        ia = 0
        for (i in 1:I) {
            for (j in 1:sec[i]) {
                ehattmp = ehat[(ic + 1):(ic + mat[i, j])]
                ic = ic + mat[i, j]
                uhatij[i, j] = huber(ehattmp)$mu
            }
            uhati[i] = huber(uhatij[i, 1:sec[i]])$mu
            for (j in 1:sec[i]) {
                reij[i, j] = uhatij[i, j] - uhati[i]
            }
            for (j in 1:sec[i]) {
                ehat2[(id + 1):(id + mat[i, j])] = ehat[(id + 
                  1):(id + mat[i, j])] - reij[i, j]
                id = id + mat[i, j]
            }
            rei[i] = huber(ehat2[(ia + 1):(ia + ni[i])])$mu
            ia = ia + ni[i]
        }
        frei = rei - huber(rei)$mu
        ib = 0
        for (i in 1:I) {
            for (j in 1:sec[i]) {
                ib = ib + 1
                rew[ib] = reij[i, j]
            }
        }
        frew = rew - huber(rew)$mu
    }
    if (location == 4) {
        matre = matrix(c(0, 0), ncol = 2)
        uhatij = 0 * mat
        ni = apply(mat, 1, sum)
        reij = 0 * mat
        rei = rep(0, I)
        rew = rep(0, sum(sec))
        frei = rep(0, I)
        uhati = rep(0, I)
        ehat2 = ehat
        ehat3 = 0 * ehat
        ic = 0
        id = 0
        ia = 0
        for (i in 1:I) {
            for (j in 1:sec[i]) {
                ehattmp = ehat[(ic + 1):(ic + mat[i, j])]
                ic = ic + mat[i, j]
                uhatij[i, j] = mean(ehattmp)
            }
            uhati[i] = mean(uhatij[i, 1:sec[i]])
            for (j in 1:sec[i]) {
                reij[i, j] = uhatij[i, j] - uhati[i]
            }
            for (j in 1:sec[i]) {
                ehat2[(id + 1):(id + mat[i, j])] = ehat[(id + 
                  1):(id + mat[i, j])] - reij[i, j]
                id = id + mat[i, j]
            }
            rei[i] = mean(ehat2[(ia + 1):(ia + ni[i])])
            ia = ia + ni[i]
        }
        frei = rei - mean(rei)
        ib = 0
        for (i in 1:I) {
            for (j in 1:sec[i]) {
                ib = ib + 1
                rew[ib] = reij[i, j]
            }
        }
        frew = rew - mean(rew)
    }
    if (scale == 1) {
        siga2 = mad(frei)^2
        sigw2 = mad(frew)^2
        ic = 0
        id = 0
        for (i in 1:I) {
            re1 = frei[i]
            for (j in 1:sec[i]) {
                id = id + 1
                re2 = frew[id]
                for (k in 1:mat[i, j]) {
                  ic = ic + 1
                  ehat3[ic] = ehat[ic] - re1 - re2
                }
            }
        }
        sigmae2 = mad(ehat3)^2
    }
    if (scale == 2) {
        siga2 = dispvar(frei)^2
        sigw2 = dispvar(frew)^2
        ic = 0
        id = 0
        for (i in 1:I) {
            re1 = frei[i]
            for (j in 1:sec[i]) {
                id = id + 1
                re2 = frew[id]
                for (k in 1:mat[i, j]) {
                  ic = ic + 1
                  ehat3[ic] = ehat[ic] - re1 - re2
                }
            }
        }
        sigmae2 = dispvar(ehat3)^2
    }
    if (scale == 3) {
        siga2 = onesampwil(frei, maktable = F, plotb = F)$tau^2
        sigw2 = onesampwil(frew, maktable = F, plotb = F)$tau^2
        ic = 0
        id = 0
        for (i in 1:I) {
            re1 = frei[i]
            for (j in 1:sec[i]) {
                id = id + 1
                re2 = frew[id]
                for (k in 1:mat[i, j]) {
                  ic = ic + 1
                  ehat3[ic] = ehat[ic] - re1 - re2
                }
            }
        }
        sigmae2 = onesampwil(ehat3, maktable = F, plotb = F)$tau^2
    }
    if (scale == 4) {
        siga2 = sd(frei)^2
        sigw2 = sd(frew)^2
        ic = 0
        id = 0
        for (i in 1:I) {
            re1 = frei[i]
            for (j in 1:sec[i]) {
                id = id + 1
                re2 = frew[id]
                for (k in 1:mat[i, j]) {
                  ic = ic + 1
                  ehat3[ic] = ehat[ic] - re1 - re2
                }
            }
        }
        sigmae2 = sd(ehat3)^2
    }
    if (scale == 5) {
        siga2 = huber(frei)$s^2
        sigw2 = huber(frew)$s^2
        ic = 0
        id = 0
        for (i in 1:I) {
            re1 = frei[i]
            for (j in 1:sec[i]) {
                id = id + 1
                re2 = frew[id]
                for (k in 1:mat[i, j]) {
                  ic = ic + 1
                  ehat3[ic] = ehat[ic] - re1 - re2
                }
            }
        }
        sigmae2 = huber(ehat3)$s^2
    }
    sigmae2 = sigmae2
    list(frei = frei, frew = frew, free = ehat3, siga2 = siga2, 
        sigw2 = sigw2, sigmae2 = sigmae2)
}
