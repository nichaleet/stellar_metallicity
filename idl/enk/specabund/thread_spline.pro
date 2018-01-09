function thread_spline, lambda, spec, ivar, mask=mask, bkspace=bkspace
    useinterp = 1

    n = n_elements(lambda)
    if ~keyword_set(mask) then mask0 = lonarr(n)+1 else mask0 = mask
    if ~keyword_set(bkspace) then bkspace = 150

    mask0[0:3] = 0
    mask0[n-4:n-1] = 0
    if n eq 8192 then mask0[4092:4100] = 0

    won = where(mask0 eq 1, con)
    bkpt = slatec_splinefit(lambda[won], spec[won], coeff, invvar=ivar[won], bkspace=bkspace, upper=1, lower=1, /silent)
    if (size(bkpt))[0] eq 0 then begin
        message, 'THREAD_SPLINE failed.', /info
        return, spec*0d + 1d
    endif
    cont = slatec_bvalu(lambda, bkpt, coeff)
        
    return, cont
end
