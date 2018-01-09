function get_sps_obs_sdss, xin, a
    common sps_spec, sps, spsz, spsage
    common get_sps, dlam, dataivar, datalam, wonfit, contmask, normalize,rest

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

    ;smooth to data wavelengths
    spsspec = smooth_gauss_wrapper(lambda*(redshift+1.), spsspec, datalam, sqrt(dlam^2+(vdisp/clight*datalam)^2))
    lambda  = datalam ;datalam is science.lambda
  
    if normalize eq 1 or rest eq 1 then begin
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;fit continuum to synthetic spectra
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       spsspec_original = spsspec
       won = where(contmask eq 1, con)
       if con lt 25 then message, 'Not enough pixels.'
       invvar=dataivar[won]

       bkpt = slatec_splinefit(lambda[won], spsspec[won], coeff, bkspace=165,invvar=invvar,upper=3, lower=3, /silent) ;SDSS   
       if bkpt[0] eq -1 then message, 'Could not fit a spline to spsspec.'
       cont = slatec_bvalu(lambda, bkpt, coeff)
       spsspec /= cont
    endif

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    if rest eq 0 then begin
       ;;only return the unmasked range (same elements as xin)
       spsspec = spsspec[wonfit]
       if total(xin-lambda(wonfit)) ne 0 then stop
       return, spsspec
    endif

    if rest eq 1 then begin ; This is to replace the get_sps_rest.pro
       if total(xin-lambda) ne 0 then stop
       return,[[spsspec],[spsspec_original],[cont]]
    endif
end
