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

    readcol, getenv('UCI')+'sps_fit/lines.txt', linestart, lineend, linetype, format='D,D,A,X', /silent, comment='#'
    contmask = bytarr(n_elements(lambda))+1
    for i=0,n_elements(linestart)-1 do begin
        w = where(lambda ge linestart[i] and lambda le lineend[i], c)
        if c gt 0 then contmask[w] = 0
    endfor
    won = where(contmask eq 1, con)
    if con lt 25 then message, 'Not enough pixels.'
    bkpt = slatec_splinefit(lambda[won], spsspec[won], coeff, bkspace=330, upper=50, lower=0.01, /silent)
    if bkpt[0] eq -1 then message, 'Could not fit a spline to spsspec.'
    cont = slatec_bvalu(lambda, bkpt, coeff)
    spsspec /= cont

    clight = 299792.458
    spsspec = smooth_gauss_wrapper(lambda, spsspec, lambda, vdisp/clight/2.35*lambda)
    spsspec = smooth_gauss_wrapper(lambda, spsspec, xin, dlam)
    return, spsspec
end
