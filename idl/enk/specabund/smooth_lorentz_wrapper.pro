function smooth_lorentz_wrapper, lambda1, spec1, lambda2, dlam, ivar1=ivar1, ivar2=ivar2
    n1 = n_elements(lambda1)
    n2 = n_elements(lambda2)
    f = long(findex(lambda1, lambda2))
    spec2 = dblarr(n2)
    if n_elements(dlam) eq 1 then dlam = replicate(dlam, n2)
    dlambda1 = -1*ts_diff(lambda1, 1)
    dlambda1 = dlambda1[where(dlambda1 gt 0)]
    halfwindow = ceil(1.1*5.*max(dlam)/min(dlambda1))

    soname = getenv('ENK_IDL')+'specabund/smooth_lorentz'+strtrim(!VERSION.MEMORY_BITS, 2)+'.so'
    
    ;now = systime(1)
    if keyword_set(ivar1) then begin
        ivar2 = dblarr(n2)
        rr = call_external(soname, 'smooth_lorentz', double(lambda1), double(spec1), double(lambda2), double(spec2), double(dlam), double(ivar1), double(ivar2), f, n1, n2, halfwindow, /i_value, value=[0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1])
    endif else begin
        rr = call_external(soname, 'smooth_lorentz_noivar', double(lambda1), double(spec1), double(lambda2), double(spec2), double(dlam), f, n1, n2, halfwindow, /i_value, value=[0, 0, 0, 0, 0, 0, 1, 1, 1])
    endelse
    ;print, systime(1)-now

   return, spec2
end
