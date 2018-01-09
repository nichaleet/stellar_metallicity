function get_sps_obs_choi, xinn, a, bkspace
    common sps_spec, sps, spsz, spsage
    common get_sps, dlam, dataivar, datalam, wonfit, npoly, contmask, normalize,rest
    z = a[0]
    age = a[1]
    vdisp = a[2]
    redshift = a[3]
    
    ;spsspec = dblarr(n_elements(xin)) - 99999.
    ;if z lt min(spsz) or z gt max(spsz) then return, spsspec
    ;if age lt min(spsage) or age gt max(spsage) then return, spsspec

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
    lambda  = lambda*(redshift+1d)       ;change lambda to observed lambda

    ;smooth to data wavelengths
    spsspec = smooth_gauss_wrapper(lambda, spsspec, lambda, vdisp/clight/2.35*lambda)
    spsspec = smooth_gauss_wrapper(lambda, spsspec, datalam, dlam)

    lambda  = datalam
  
    ;fit continuum to synthetic spectra
    spline = 0
    poly = 0
    simpoly = 1

    readcol, '/scr2/nichal/workspace2/sps_fit/lines.txt', linestart, lineend, linetype, format='D,D,A,X', /silent, comment='#'
    contmask = bytarr(n_elements(lambda))+1
    for i=0,n_elements(linestart)-1 do begin
       if linetype[i] eq 'a' then begin
          w = where(lambda/(redshift+1d) ge linestart[i] and lambda/(redshift+1d) le lineend[i], c)
        ;  if c gt 0 then contmask[w] = 0
       endif
     endfor

    won = where(contmask eq 1, con)
    if con lt 25 then message, 'Not enough pixels.'
    invvar=dataivar[won]/(median(spsspec[won]))^2

    case 1 of
       spline: begin
          bkpt = slatec_splinefit(lambda[won], spsspec[won], coeff, bkspace=200,invvar=invvar,upper=5, lower=1.5, /silent,/everyn,mask=mask) ;SDSS
   
          if bkpt[0] eq -1 then message, 'Could not fit a spline to spsspec.'
          cont = slatec_bvalu(lambda, bkpt, coeff)
       end
       poly: begin
          degree = 6
          norm = median(spsspec[won])
          a = [norm, replicate(0.0, degree-1)]
          p = lmfit(lambda[won],spsspec[won], a, measure_errors=(invvar)^(-0.5), /double, function_name='legendre_poly')
          cont = legendre_poly(lambda, a, /noderiv)
       end
       simpoly:begin
          degree = npoly
          divfac = median(spsspec[won])
          p=poly_fit(lambda[won],spsspec[won]-divfac,degree)
          cont = poly(lambda,p)+divfac
       end
    endcase
    spsspec /= cont
    
    ;only return the unmasked range (same elements as xin)
    spsspec = spsspec[wonfit]
    if total(xinn-lambda(wonfit)) ne 0 then stop

   ; plot, lambda,spsspec,yrange=[0.8,1.2],xrange=[3000,7000]
    ;spsspec = interpol(spsspec,lambda,xin,/lsquadratic)
    ;oplot, lambda,spsspec,color=255
    ;wait,0.5
    return, spsspec
end
