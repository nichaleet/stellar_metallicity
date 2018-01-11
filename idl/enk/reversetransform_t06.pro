function reversetransform_t06, sdss, sdsserr, error=error
    ;ugriz -> u'g'r'i'z'
    ;Tucker et al. (2006)

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
    sdss = dblarr(5)
    error = dblarr(5)
    sdss[0] = u
    sdss[1] = 0.943396*g + 0.054611*r + 0.001993*i + 0.041887
    sdss[2] = r - 0.035211*(r-i-0.21)
    sdss[3] = i - 0.041248*(r-i-0.21)
    sdss[4] = 0.970874*z + 0.030327*i - 0.001201*r - 0.008737
    error[0] = uerr
    error[1] = sqrt((0.943396*gerr)^2. + (0.054611*rerr)^2. + (0.001993*ierr)^2.)
    error[2] = sqrt((0.964789*rerr)^2. + (0.035211*ierr)^2.)
    error[3] = sqrt((0.041248*rerr)^2. + (1.041248*ierr)^2.)
    error[4] = sqrt((0.970874*zerr)^2. + (0.030327*ierr)^2. + (0.001201*rerr)^2.)
    return, sdss
end
