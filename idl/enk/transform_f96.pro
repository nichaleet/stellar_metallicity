function transform_f96, sdss, sdsserr, error=error
    ;u'g'r'i'z' -> UBVRI
    ;Fukugita et al. 1996
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
