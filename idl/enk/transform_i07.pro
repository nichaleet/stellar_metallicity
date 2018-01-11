function transform_i07, sdss, sdsserr, error=error
    ;ugriz -> UBVRI
    ;Ivezic et al. (2007)
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
    ubvri = dblarr(5)
    error = dblarr(5)
    ubvri[0] = -999d
    ubvri[1] = -999d
    ubvri[2] = g - 0.56*(g-r+0.23)/1.05 + 0.12
    ubvri[3] = (r - ubvri[2] - 0.13)/0.84 + ubvri[2]
    ubvri[4] = ubvri[3] - (r-i+0.23)/0.98
    error[0] = -999d
    error[1] = -999d
    error[2] = 0d
    error[3] = -999d
    error[4] = 0d
    return, ubvri
end
