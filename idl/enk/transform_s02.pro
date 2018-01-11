function transform_s02, sdss, sdsserr, error=error
    ;u'g'r'i'z' -> UBVRI
    ;Smith et al. (2002) AJ, 123, 2121S
    sdssprime = reversetransform_t06(sdss, sdsserr, error=sdssprimeerr)
    u = sdssprime[0]
    g = sdssprime[1]
    r = sdssprime[2]
    i = sdssprime[3]
    z = sdssprime[4]
    uerr = sdssprimeerr[0]
    gerr = sdssprimeerr[1]
    rerr = sdssprimeerr[2]
    ierr = sdssprimeerr[3]
    zerr = sdssprimeerr[4]
    ubvri = dblarr(5)
    error = dblarr(5)
    ubvri[0] = 0.75*u + 0.72*g - 0.47*r - 0.66
    ubvri[1] = g + 0.47*(g-r) + 0.17
    ubvri[2] = g - 0.55*(g-r) - 0.03
    ubvri[3] = -0.14*g + 1.14*r - 0.14
    ubvri[4] = r-i lt 0.95 ? -0.14*g + 0.14*r + 1.00*i - 0.35 : -0.14*g + 0.44*r +0.70*i - 0.63
    error[0] = sqrt((0.75*uerr)^2. + (0.72*gerr)^2. + (0.47*rerr)^2.)
    error[1] = sqrt((1.47*gerr)^2. + (0.47*rerr)^2.)
    error[2] = sqrt((0.45*gerr)^2. + (0.55*rerr)^2.)
    error[3] = sqrt((0.14*gerr)^2. + (1.14*rerr)^2.)
    error[4] = r-i lt 0.95 ? sqrt((0.14*gerr)^2. + (0.14*rerr)^2. + (1.00*ierr)^2.) : sqrt((0.14*gerr)^2. + (0.44*rerr)^2. + (0.70*ierr)^2.)
    return, ubvri
end
