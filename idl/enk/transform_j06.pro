function transform_j06, sdss, sdsserr, error=error
    ;ugriz -> UBVRI
    ;Jordi et al. (2006)
    u = sdss[0]
    g = sdss[1]
    r = sdss[2]
    i = sdss[3]
    z = sdss[4]
    ugood = u gt 0.0 and u lt 35.0
    ggood = g gt 0.0 and g lt 35.0
    rgood = r gt 0.0 and r lt 35.0
    igood = i gt 0.0 and i lt 35.0
    zgood = z gt 0.0 and z lt 35.0
    uerr = sdsserr[0]
    gerr = sdsserr[1]
    rerr = sdsserr[2]
    ierr = sdsserr[3]
    zerr = sdsserr[4]
    ubvri = dblarr(5)
    error = dblarr(5)

    ;global
    a1 = 0.630  &  a1err = 0.002  &  b1 = -0.124  &  b1err = 0.002
    a2 = 1.007  &  a2err = 0.005  &  b2 = -0.236  &  b2err = 0.003
    a5 = 0.750  &  a5err = 0.050  &  b5 = 0.770   &  b5err = 0.070 & c5 = 0.720 & c5err = 0.040
    a7 = 1.646  &  a7err = 0.008  &  b7 = -0.139  &  b7err = 0.004
    a8 = 0.247  &  a8err = 0.003  &  b8 = 0.329   &  b8err = 0.002
    VR = (g-r-b7)/a7
    if VR le 0.93 then begin
        a4 = 0.267  &  a4err = 0.005  &  b4 = 0.088   &  b4err = 0.003
    endif else begin
        a4 = 0.77   &  a4err = 0.04   &  b4 = -0.37   &  b4err = 0.04
    endelse

    ;;Population I
    ;a1 = 0.634  &  a1err = 0.002  &  b1 = -0.127  &  b1err = 0.002
    ;a2 = 0.988  &  a2err = 0.006  &  b2 = -0.221  &  b2err = 0.004
    ;a7 = 1.599  &  a7err = 0.009  &  b7 = -0.106  &  b7err = 0.006
    ;a8 = 0.251  &  a8err = 0.003  &  b8 = 0.335   &  b8err = 0.002
    ;VR = (g-r-b7)/a7
    ;if VR le 0.93 then begin
    ;    a4 = 0.275  &  a4err = 0.006  &  b4 = 0.086   &  b4err = 0.004
    ;endif else begin
    ;    a4 = 0.71   &  a4err = 0.05   &  b4 = -0.31   &  b4err = 0.05
    ;endelse

    ;;metal-poor Population II
    ;a1 = 0.596  &  a1err = 0.009  &  b1 = -0.148  &  b1err = 0.007
    ;a2 = 1.06   &  a2err = 0.02   &  b2 = -0.30   &  b2err = 0.01
    ;a4 = 0.34   &  a4err = 0.02   &  b4 = 0.015   &  b4err = 0.008
    ;a7 = 1.72   &  a7err = 0.02   &  b7 = -0.198  &  b7err = 0.007
    ;a8 = 0.21   &  a8err = 0.02   &  b8 = 0.34    &  b8err = 0.01

    ubvri[0] = ugood and ggood and rgood ? (u-g - b5*(g-(r - (a4-1d)*(g-r-b7)/a7 - b4)-b1)/a1 - c5)/a5 + (g + (a1-1d)*(r - (a4-1d)*(g-r-b7)/a7 - b4) - b1) / a1 : -999d
    ubvri[1] = ggood and rgood ? (g + (a1-1d)*(r - (a4-1d)*(g-r-b7)/a7 - b4) - b1) / a1 : -999d
    ubvri[2] = ggood and rgood ? r - (a4-1d)*(g-r-b7)/a7 - b4 : -999d
    ubvri[3] = rgood and igood ? i - (a8-1d)*(r-i-b2)/a2 - b8 : -999d
    ubvri[4] = rgood and igood ? i - a8*(r-i-b2)/a2 - b8 : -999d
    error[0] = ugood and ggood and rgood ? sqrt((uerr/a5)^2. + $
                    (gerr*(1d/a1 - (a1-1d)*(a4-1d)/(a1*a7) - 1d/a5 - b5/(a1*a5) - b5*(a4-1d)/(a1*a5*a7)))^2. + $
                    (rerr*(b5/(a1*a5) + b5*(a4-1d)/(a1*a5*a7) + (a4-1d)/a1))^2. + $
                    (a1err*(b5*(g-r-(a4-1d)*(g-r-b7)/a7+b4-b1)/(a5*a1^2.) + (r-(a4-1d)*(g-r-b7)/a7-b4)/a1 - (g+(a1-1d)*(r-(a4-1d)*(g-r-b7)/a7-b4)-b1)/a1^2.))^2. + $
                    (a4err*(-1d*b5*(g-r-b7)/(a1*a5*a7) - (a1-1d)*(g-r-b7)/(a1*a7)))^2. + $
                    (a5err*(-1d*(u-g-b5*(g-r+(a4-1d)*(g-r-b7)/a7-b1)/a1-c5)/a5^2.))^2. + $
                    (a7err*(-1d*b5*(a4-1d)*(g-r-b7)/(a1*a5*a7^2.) - (a1-1d)*(a4-1d)*(g-r-b7)/(a1*a7^2.)))^2. + $
                    (b1err*(b5/(a1*a5) + (1d/a1)))^2. + $
                    (b4err*(-1d*b5/(a1*a5) - (a1-1d)/a1))^2. + $
                    (b5err*(-1d*(g-r+(a4-1d)*(g-r-b7)/a7+b4-b1)/(a1*a5)))^2. + $
                    (b7err*(b5*(a4-1d)/(a1*a5*a7) + (a1-1d)*(a4-1d)/(a1*a7)))^2. + $
                    (c5err*(-1d/a5))^2.) : -999d
    error[1] = ggood and rgood ? sqrt((gerr*(1d/a1) - (a1-1d)*(a4-1d)/(a1*a7))^2. + $
                    (rerr*((a1-1d)/a1 + (a4-1d)/(a1*a7)))^2. + $
                    (a1err*((r-(a4-1d)*(g-r-b7)/a7-b4)/a1 - (g+(a1-1d)*(r-(a4-1)*(g-r-b7)/a7-b4)-b1)/(a1^2.)))^2. + $
                    (a4err*(-1d*(a1-1d)*(g-r-b7)/(a1*a7)))^2. + $
                    (a7err*(a1-1d)*(a4-1d)*(g-r-b7)/(a1*a7^2.))^2. + $
                    (b1err/a1)^2. + $
                    (b4err*(-1d)*(a1-1d)/a1)^2. + $
                    (b7err*(a1-1d)*(a4-1d)/(a1*a7))^2.) : -999d
    error[2] = ggood and rgood ? sqrt((gerr*(-1d)*(a4-1d)/a7)^2. + $
                    (rerr*(1d + (a4-1d)/a7))^2. + $
                    (a4err*(-1d)*(g-r-b7)/a7)^2. + $
                    (a7err*(a4-1d)*(g-r-b7)/(a7^2.))^2. + $
                    (b4err*(-1d))^2. + $
                    (b7err*(a4-1d)/a7)^2.) : -999d
    error[3] = rgood and igood ? sqrt((rerr*(-1d)*(a8-1d)/a2)^2. + $
                    (ierr*(1d + (a8-1d)/a2))^2. + $
                    (a2err*(a8-1d)*(r-i-b2)/(a2^2.))^2. + $
                    (a8err*(-1d)*(r-i-b2)/a2)^2. + $
                    (b2err*(a8-1d)/a2)^2. + $
                    (b8err*(-1d))^2.) : -999d
    error[4] = rgood and igood ? sqrt((rerr*(-1d)*a8/a2)^2. + $
                    (ierr*(1d + a8/a2))^2. + $
                    (a2err*a8*(r-i-b2)/(a2^2.))^2. + $
                    (a8err*(-1d)*(r-i-b2)/a2)^2. + $
                    (b2err*a8/a2)^2. + $
                    (b8err*(-1d))^2.) : -999d
    return, ubvri
end
