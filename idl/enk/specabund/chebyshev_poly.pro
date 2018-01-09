function chebyshev_poly, x, a, kind=kind, noderiv=noderiv
    if ~keyword_set(kind) then kind = 1 else begin
        if kind ne 1 and kind ne 2 then message, 'Chebyshev polynomials must be of the first or second kind.'
    endelse    
    n = n_elements(a)
    Ti2 = x*0d0 + 1d0
    y = a[0]*Ti2
    if ~keyword_set(noderiv) then begin
        dyda = dblarr(n, n_elements(x))
        dyda[0,*] = Ti2
    endif
    if n eq 1 then return, keyword_set(noderiv) ? y : [transpose(reform(y)), dyda]
    Ti1 = double(kind)*x
    y += a[1]*Ti1
    if ~keyword_set(noderiv) then dyda[1,*] = Ti1
    if n eq 2 then return, keyword_set(noderiv) ? y : [transpose(reform(y)), dyda]
    for i=2,n-1 do begin
        Ti = 2d0*x*Ti1 - Ti2
        y += a[i]*Ti
        if ~keyword_set(noderiv) then dyda[i,*] = Ti
        Ti2 = Ti1
        Ti1 = Ti
    endfor
    return, keyword_set(noderiv) ? y : [transpose(reform(y)), dyda]
end
