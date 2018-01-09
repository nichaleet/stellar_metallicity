pro specres,fit,str
    npix = n_elements(str.lambda)
    nsky = 200
    lambda = str.lambda
    str=create_struct(str,'dlam',dblarr(npix))
    str=create_struct(str,'skyfit',[[dblarr(nsky)], [dblarr(nsky)], [dblarr(nsky)]])
    str=create_struct(str,'skylinemask',lonarr(nsky))
    str=create_struct(str,'goodsky',0)

    n = (size(fit))[1]
    if n gt 200 then message, 'Too many sky lines!'
    if n gt 0 then str.skyfit[0:n-1,*] = fit
    if n le 3 then begin
       str.dlam = replicate(1.37/2.35, n_elements(lambda))
       str.skylinemask = lonarr(n_elements(str.skylinemask))-1
       str.goodsky = 0
       message, 'Unstable sky line fit.  Using FWHM = 1.37 A', /info
       return
    endif

    w = where(2.35*fit[*,1] gt 0.8 and 2.35*fit[*,1] lt 7.0, cwprev)
    if cwprev lt 3 then begin
       str.dlam = replicate(1.37/2.35, n_elements(lambda))
       str.skylinemask = lonarr(n_elements(str.skylinemask))-1
       str.goodsky = 0
       message, 'Unusuable arc lines.  Using FWHM = 1.37 A', /info
       return
    endif

    ;quadratic fit
    qf = poly_fit(fit[w,0]/1000.0 - 7.8, 2.35*fit[w,1], 2, measure_errors=2.35*fit[w,2], chisq=chisq, /double, yfit=yfit)
    ;remove outliers and refit 4 times
    for j=0,4 do begin
       wnow = where(abs(2.35*fit[w,1] - yfit) lt 2.*2.35*fit[w,2], cw) ;good sigmas
       if cw eq cwprev then break
       cwprev = cw
       if cw lt 3 then begin
          str.goodsky = 0
          message, 'The spectral resolution fit is very poor.', /info
          break
       endif
       w = w[wnow]
       qf = poly_fit(fit[w,0]/1000.0 - 7.8, 2.35*fit[w,1], 2, measure_errors=2.35*fit[w,2], chisq=chisq, /double, yfit=yfit)
    endfor

    n = (size(fit))[1]
    str.skylinemask = 0
    str.skylinemask[w] = 1

    if n lt 200 then str.skylinemask[n:n_elements(str.skylinemask)-1] = -1

    l = lambda / 1000. - 7.8
    dlam = poly(l, qf)
    dlam /= 2.35
    str.dlam = dlam   ;dlam is not FWHM. It is sigma!!
end
