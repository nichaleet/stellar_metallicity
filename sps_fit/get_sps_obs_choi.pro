function get_sps_obs_choi, xin, a, bkspace
    common sps_spec, sps, spsz, spsage
    common get_sps, dlam, dataivar, datalam, wonfit, npoly, contmask, normalize,rest
    common get_coelho,coelho_normalize

    z = a[0]
    age = a[1]
    vdisp = a[2]
    redshift = a[3]
    if n_Elements(a) eq 5 then begin
       yesalpha=1
       afe = a[4]
    endif else yesalpha=0
    
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
    lambda  = datalam           ;datalam is science.lambda
      

    ;add response function to [a/Fe]
    if yesalpha then begin
       response_fn_str = coelho_response_fn(z,age,[afe],vdisp,redshift)
       response_fn = response_fn_str.dspec
       mgbreg = where(lambda/(1.+redshift) gt 5160. and lambda/(1.+redshift) lt 5192.,cmgbreg)
       if cmgbreg lt 5 then message,'Mgb Region is too small'
       spsspec(mgbreg) = spsspec(mgbreg)*response_fn(mgbreg)
       ;spsspec = spsspec*response_fn
    endif

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;fit continuum to synthetic spectra
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    if normalize eq 1 or rest eq 1 then begin
       spline = 1
       poly = 0
       simpoly = 0

       spsspec_original = spsspec
       won = where(contmask eq 1, con)
       if con lt 25 then message, 'Not enough pixels.'
       invvar=dataivar[won]

       case 1 of
          spline: begin
             bkpt = slatec_splinefit(lambda[won], spsspec[won], coeff, bkspace=165,invvar=invvar,upper=3, lower=3, /silent,/everyn,mask=mask) ;SDSS
             
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
             p=poly_fit(lambda[won],spsspec[won],degree)
             cont = poly(lambda,p)
          end
       endcase
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
