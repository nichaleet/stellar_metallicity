function transform_cg08, sdss, sdsserr, error=error
    ;ugriz -> UBVRI
    ;Chonis & Gaskill 2008
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
    ubvri[0] = u - 0.854
    ubvri[1] = g + 0.327*(g-r) + 0.216
    ubvri[2] = g - 0.587*(g-r) - 0.011
    ubvri[3] = r - 0.272*(r-i) - 0.159
    ubvri[4] = i - 0.337*(r-i) - 0.370
    error[0] = sqrt(uerr^2. + 0.007^2.)
    error[1] = sqrt((1.327*gerr)^2. + (0.327*rerr)^2. + (0.047*(g-r))^2. + 0.027^2.)
    error[2] = sqrt((0.413*gerr)^2. + (0.587*rerr)^2. + (0.022*(g-r))^2. + 0.013^2.)
    error[3] = sqrt((0.728*rerr)^2. + (0.272*ierr)^2. + (0.092*(r-i))^2. + 0.022^2.)
    error[4] = sqrt((0.663*ierr)^2. + (0.337*rerr)^2. + (0.191*(r-i))^2. + 0.041^2.)
    return, ubvri
end
