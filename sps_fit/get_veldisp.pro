function get_veldisp, xin, a
    common get_sps, dlam, dataivar, datalam, wonfit, contmask, normalize

    redshift = a[0]
    vdisp = a[1]

    spsstruct = sps_interp(z, age)
    lambda = spsstruct.lambda
    spsspec = spsstruct.spec

    w = where(lambda gt 0 and lambda lt 100000, c)
    if c lt 25 then message, 'Not enough pixels.'
    lambda = lambda[w]
    spsspec = spsspec[w]
    clight = 299792.458

    spsspec = spsspec*clight/lambda^2    ;change fnu(Lsun/Hz) to flambda
    spsspec = spsspec/median(spsspec)    ;normalize to around 1

    ;smooth to data wavelengths
    spsspec = smooth_gauss_wrapper(lambda, spsspec, lambda, vdisp/clight/2.35*lambda)
    spsspec = smooth_gauss_wrapper(lambda*(redshift+1.), spsspec, datalam, dlam)

    lambda  = datalam ;datalam is science.lambda
  
    if normalize eq 1 then begin
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;fit continuum to synthetic spectra
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       n = n_elements(lambda)
       nhalf = round(double(n)/2.)
       if n lt 8193 then begin
          wwhole = lindgen(nhalf)
          ccd2 = 2
       endif else begin
          wwhole = lindgen(n)
          ccd2 = 1
       endelse
       for ccd=1,ccd2 do begin
          won = where(contmask[wwhole] eq 1, complement=woff, con) + wwhole[0]
          woff += wwhole[0]   
          if con lt 25 then message, 'Not enough pixels.'
          invvar = dataivar[won]/(median(spsspec[won]))^2
          bkpt = slatec_splinefit(lambda[won], spsspec[won], coeff, invvar=invvar, bkspace=150, upper=3, lower=3, /silent) ;DEIMOS
          if bkpt[0] eq -1 then message, 'Could not fit a spline to spsspec.'
          cont = slatec_bvalu(lambda[wwhole], bkpt, coeff)
          
          if ccd eq 1 then contb = cont
          if ccd eq 2 then contr = cont
          wwhole += nhalf
       endfor
       if n lt 8193 then cont = [contb, contr] else cont = contb
       
       spsspec /= cont
    endif
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;only return the unmasked range (same elements as xin)
    spsspec = spsspec[wonfit]
    if total(xin-lambda(wonfit)) ne 0 then stop

    return, spsspec
end
