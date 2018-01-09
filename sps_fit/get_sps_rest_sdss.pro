function get_sps_rest_sdss, xin, a
    common sps_spec, sps, spsz, spsage
    common get_sps, dlam, dataivar, datalam, wonfit

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

    spsspec = smooth_gauss_wrapper(lambda, spsspec, lambda, vdisp/clight/2.35*lambda)
    spsspec = smooth_gauss_wrapper(lambda, spsspec, datalam, dlam)
    lambda  = datalam
 

    spsspec_original = spsspec
    readcol, '/scr2/nichal/workspace2/sps_fit/lines.txt', linestart, lineend, linetype, format='D,D,A,X', /silent, comment='#'
    contmask = bytarr(n_elements(lambda))+1
    for i=0,n_elements(linestart)-1 do begin
       w = where(lambda ge linestart[i] and lambda le lineend[i], c)
       if c gt 0 then contmask[w] = 0
    endfor

    won = where(contmask eq 1, con)
    ;print,'ncon = ',con,' ',n_elements(contmask)
    if con lt 25 then message, 'Not enough pixels.'
    ;bkpt = slatec_splinefit(lambda[won], spsspec[won], coeff, bkspace=330, upper=5, lower=1, /silent,/everyn);DEIMOS
    invvar=dataivar[won]/(median(spsspec[won]))^2
  
    bkpt = slatec_splinefit(lambda[won], spsspec[won], coeff, bkspace=165,invvar=invvar, upper=5, lower=1.5, /silent,/everyn,mask=mask) ;SDSS 
    if bkpt[0] eq -1 then message, 'Could not fit a spline to spsspec.'
    cont = slatec_bvalu(lambda, bkpt, coeff)    

    spsspec /= cont

    if total(xin-lambda(wonfit)) ne 0 then stop

    return, [[spsspec[wonfit]],[spsspec_original[wonfit]],[cont[wonfit]]]
end
