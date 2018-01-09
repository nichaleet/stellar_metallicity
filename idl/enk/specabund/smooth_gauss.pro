function smooth_gauss, lambda1, spec1, lambda2, dlam, ivar1=ivar1, ivar2=ivar2, maxsigma=maxsigma
    n1 = n_elements(lambda1)
    n2 = n_elements(lambda2)
    f = long(findex(lambda1, lambda2))
    spec2 = dblarr(n2)
    ivar2 = dblarr(n2)
    if n_elements(dlam) eq 1 then dlam = replicate(dlam, n2)
    if ~keyword_set(ivar1) then ivar1 = replicate(1.0, n1)
    dlambda1 = -1*ts_diff(lambda1, 1)
    dlambda1 = dlambda1[where(dlambda1 gt 0)]
    if ~keyword_set(maxsigma) then maxsigma = 4.
    halfwindow = ceil(1.1*maxsigma*max(dlam)/min(dlambda1))
    ;acw = lonarr(n2)
    ;wheretime = 0d
    ;exptime = 0d
    ;totaltime = 0d
    ;now = systime(1)
    for i=0L,n2-1 do begin
        ;now = systime(1)
        low = (f[i]-halfwindow) > 0
        high = ((f[i]+halfwindow) < (n1-1)) > 0
        if low lt n1 and low lt high then begin
            w = where(abs(lambda1[low:high] - lambda2[i]) lt dlam[i]*maxsigma, cw) + low
            ;acw[i] = cw
            ;wheretime += -1*now + (now = systime(1))
            if cw gt 0 then begin
                gauss = exp(-1.0*(lambda1[w] - lambda2[i])^2 / (2.0*dlam[i]^2))
                ;print, total(gauss)/(dlam[i]*sqrt(2.*!DPI))*mean(dlambda1)
                ;exptime += -1*now + (now = systime(1))
                temp = ivar1[w]*gauss
                temp2 = total(temp)
                spec2[i] = total(spec1[w]*temp) / temp2
                ivar2[i] = temp2^2. / total(temp*gauss)
                ;totaltime += -1*now + (now = systime(1))
            endif
        endif
    endfor
    ;print, wheretime, exptime, totaltime, wheretime+exptime+totaltime
    ;print, systime(1)-now
    return, spec2
end
