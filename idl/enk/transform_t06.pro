function transform_t06, ugrizprime, ugrizprimeerr, error=sdsserr
    ;u'g'r'i'z' -> ugriz
    ;Tucker et al. (2006)

    u = ugrizprime[0]
    g = ugrizprime[1]
    r = ugrizprime[2]
    i = ugrizprime[3]
    z = ugrizprime[4]
    uerr = ugrizprimeerr[0]
    gerr = ugrizprimeerr[1]
    rerr = ugrizprimeerr[2]
    ierr = ugrizprimeerr[3]
    zerr = ugrizprimeerr[4]
    sdss = dblarr(5)
    sdsserr = dblarr(5)
    sdss[0] = u - 0.0556*(g-i)^3. - 0.3332*(g-i)^2. + 0.5293*(g-i) - 0.1816 ;Clem et al. (2008)
    sdss[1] = g + 0.060*(g-r-0.53)
    sdss[2] = r + 0.035*(r-i-0.21)
    sdss[3] = i + 0.041*(r-i-0.21)
    sdss[4] = z - 0.030*(i-z-0.09)
    sdsserr[0] = uerr
    sdsserr[1] = sqrt((1.060*gerr)^2. + (0.060*rerr)^2.)
    sdsserr[2] = sqrt((1.035*rerr)^2. + (0.035*ierr)^2.)
    sdsserr[3] = sqrt((1.041*rerr)^2. + (0.041*ierr)^2.)
    sdsserr[4] = sqrt((1.030*zerr)^2. + (0.030*ierr)^2.)
    return, sdss
end
