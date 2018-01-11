function transform_r09_sdss, sdss, sdsserr, error=cfhterr
    ;ugriz (SDSS) -> ugriz (CFHT)
    ;Regnault et al. (2009)

    u = sdss[0]
    g = sdss[1]
    r = sdss[2]
    i = sdss[3]
    z = sdss[4]
    uerr = sdsserr[0]
    gerr = sdsserr[1]
    rerr = sdsserr[2]
    ierr = sdsserr[3]
    zerr = sdsserr[4]
    cfht = dblarr(5)
    cfhterr = dblarr(5)
    cfht[0] = u - 0.211*(u-g) - 0.56
    cfht[1] = g - 0.155*(g-r) + 0.11
    cfht[2] = r - 0.030*(r-i) - 0.14
    cfht[3] = i - 0.102*(r-i) - 0.32
    cfht[4] = z + 0.036*(i-z) - 0.43
    cfhterr[0] = sqrt((0.789*uerr)^2. + (0.211*gerr)^2. + (0.004*(u-g))^2.)
    cfhterr[1] = sqrt((0.845*gerr)^2. + (0.155*rerr)^2. + (0.003*(g-r))^2.)
    cfhterr[2] = sqrt((0.970*rerr)^2. + (0.030*ierr)^2. + (0.004*(r-i))^2.)
    cfhterr[3] = sqrt((0.102*rerr)^2. + (0.898*ierr)^2. + (0.005*(r-i))^2.)
    cfhterr[4] = sqrt((0.036*ierr)^2. + (1.036*zerr)^2. + (0.008*(i-z))^2.)
    return, cfht
end
