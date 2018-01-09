function smooth_ang, lambda, specin, dlam
    n = n_elements(specin)
    spec = dblarr(n)
    dlam2 = dlam/2.
    for i=0L,n-1 do spec[i] = mean(specin[where(lambda gt lambda[i]-dlam2 and lambda le lambda[i]+dlam2)])
    return, spec
end
