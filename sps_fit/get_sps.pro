function get_sps, xin, a
    common sps_spec, sps, spsz, spsage
    common get_sps, dlam

    z = a[0]
    age = a[1]
    vdisp = a[2]

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
    spsspec_original = spsspec

    readcol, '/scr2/nichal/workspace2/sps_fit/lines.txt', linestart, lineend, linetype, format='D,D,A,X', /silent, comment='#'
    contmask = bytarr(n_elements(lambda))+1
    for i=0,n_elements(linestart)-1 do begin
        w = where(lambda ge linestart[i] and lambda le lineend[i], c)
        if c gt 0 then contmask[w] = 0
    endfor
    won = where(contmask eq 1, con)
    if con lt 25 then message, 'Not enough pixels.'
    bkpt = slatec_splinefit(lambda[won], spsspec[won], coeff, bkspace=330, upper=5, lower=1, /silent,/everyn)
    if bkpt[0] eq -1 then message, 'Could not fit a spline to spsspec.'
    cont = slatec_bvalu(lambda, bkpt, coeff)
    spsspec /= cont

    clight = 299792.458
    spsspec = smooth_gauss_wrapper(lambda, spsspec, lambda, vdisp/clight/2.35*lambda)
   ; plot, lambda,spsspec,yrange=[0.8,1.2],xrange=[3000,7000]
    ;spsspec = interpol(spsspec,lambda,xin,/lsquadratic)
    spsspec = smooth_gauss_wrapper(lambda, spsspec, xin, dlam)
    ;oplot, lambda,spsspec,color=255
    ;wait,0.5
    return, spsspec
end
