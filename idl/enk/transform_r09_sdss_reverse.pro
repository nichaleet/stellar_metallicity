function transform_r09_sdss_reverse, cfht, cfhterr, error=sdsserr
    ;ugriz (CFHT) -> ugriz (SDSS)
    ;Regnault et al. (2009)

    u = cfht[0]
    g = cfht[1]
    r = cfht[2]
    i = cfht[3]
    z = cfht[4]
    uerr = cfhterr[0]
    gerr = cfhterr[1]
    rerr = cfhterr[2]
    ierr = cfhterr[3]
    zerr = cfhterr[4]
    sdss = dblarr(5)
    sdsserr = dblarr(5)
    
    a0 = -0.211  &  a0err = 0.004  &  a1 = -0.56
    b0 = -0.155  &  b0err = 0.003  &  b1 = 0.11
    c0 = -0.030  &  c0err = 0.004  &  c1 = -0.14
    d0 = -0.102  &  d0err = 0.005  &  d1 = -0.32
    e0 = 0.036   &  e0err = 0.008  &  e1 = -0.43

    sdss[0] = (u + a0/(b0+1d)*(g + b0*(r + c0*(r-i-c1+d1)/(c0-d0+1d) - c1)) - a1) / (a0+1d)
    sdss[1] = (g + b0*(r + c0*(r-i-c1+d1)/(c0-d0+1d) - c1)) / (b0+1d)
    sdss[2] = r + c0*(r-i-c1+d1)/(c0-d0+1d) - c1
    sdss[3] = i - d0*(r-i-c1+d1)/(c0-d0+1d) - d1
    sdss[4] = (z - e0*(i - d0*(r-i-c1+d1)/(c0-d0+1d) - d1) - e1) / (1d - e0)

    sdsserr[0] = sqrt((uerr/(a0+1d))^2. + $
                      (gerr*a0/(b0+1d)/(a0+1d))^2. + $
                      (rerr*(a0*b0/(a0+1d)/(b0+1d) + a0*b0*c0/(a0+1d)/(b0+1d)/(c0-d0+1d)))^2. + $
                      (ierr*(-1d)*a0*b0*c0/(a0+1d)/(b0+1d)/(c0-d0+1d))^2. + $
                      (a0err*((g + b0*(r + c0*(r-i-c1+d1)/(c0-d0+1d) - c1))/(a0+1d)/(b0+1d) - (u + a0/(b0+1d)*(g + b0*(r + c0*(r-i-c1+d1)/(c0-d0+1d) - c1)) - a1) / (a0+1d)^2.))^2. + $
                      (b0err*((-1d)*a0/(b0+1d)^2.*(g + b0*(r + c0*(r-i-c1+d1)/(c0-d0+1d) - c1))/(a0+1d) + (r + c0*(r-i-c1+d1)/(c0-d0+1d) - c1)*a0/(a0+1d)/(b0+1d)))^2. + $
                      (c0err*(a0*b0/(a0+1d)/(b0+1d)*(r-i-c1+d1)/(c0-d0+1d) - a0*b0*c0/(a0+1d)/(b0+1d)*(r-i-c1+d1)/(c0-d0+1d)^2.))^2. + $
                      (d0err*(a0*b0*c0/(a0+1d)/(b0+1d)*(r-i-c1+d1)/(c0-d0+1d)^2.))^2.)
    sdsserr[1] = sqrt((gerr/(b0+1d))^2. + $
                      (rerr*(b0/(b0+1d) + b0*c0/(c0-d0+1d)/(b0+1d)))^2. + $
                      (ierr*((-1d)*b0*c0/(c0-d0+1d)/(b0+1d)))^2. + $
                      (b0err*((-1d)*(g + b0*(r + c0*(r-i-c1+d1)/(c0-d0+1d) - c1))/(b0+1d)^2. + (r + c0*(r-i-c1+d1)/(c0-d0+1d) - c1)/(b0+1d)))^2. + $
                      (c0err*(b0*(r-i-c1+d1)/(c0-d0+1d)/(b0+1d) - b0/(b0+1d)*(r-i-c1+d1)/(c0-d0+1d)^2.))^2. + $
                      (d0err*(b0*c0/(b0+1d)*(r-i-c1+d1)/(c0-d0+1d)^2.))^2.)
    sdsserr[2] = sqrt((rerr*(1d + c0/(c0-d0+1d)))^2. + $
                      (ierr*(-1d*c0/(c0-d0+1d)))^2. + $
                      (c0err*((r-i-c1+d1)/(c0-d0+1d) - c0*(r-i-c1+d1)/(c0-d0+1d)^2.))^2. + $
                      (d0err*(c0*(r-i-c1+d1)/(c0-d0+1d)^2.))^2.)
    sdsserr[3] = sqrt((rerr*(-1d*d0/(c0-d0+1d)))^2. + $
                      (ierr*(1d + d0/(c0-d0+1d)))^2. + $
                      (c0err*(d0*(r-i-c1+d1)/(c0-d0+1d)^2.))^2. + $
                      (d0err*((-1d)*(r-i-c1+d1)/(c0-d0+1d) - d0*(r-i-c1+d1)/(c0-d0+1d)^2.))^2.)
    sdsserr[4] = sqrt((rerr*(e0*d0/(c0-d0+1d)/(1d - e0)))^2. + $
                      (ierr*(e0/(e0-1d) - d0/(c0-d0+1d)/(1d - e0)))^2. + $
                      (zerr/(1d - e0))^2. + $
                      (c0err*(e0*d0/(e0-1d)*(r-i-c1+d1)/(c0-d0+1d)^2.))^2. + $
                      (d0err*(e0*(r-i-c1+d1)/(c0-d0+1d)/(1d - e0) + e0*d0/(1d - e0)*(r-i-c1+d1)/(c0-d0+1d)^2.))^2. + $
                      (e0err*((i - d0*(r-i-c1+d1)/(c0-d0+1d) - d1)/(1d - e0) - (z - e0*(i - d0*(r-i-c1+d1)/(c0-d0+1d) - d1) - e1) / (1d - e0)^2.))^2.)
    return, sdss
end
