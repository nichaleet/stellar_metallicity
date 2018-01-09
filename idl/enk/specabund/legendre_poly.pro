function legendre_poly, x, a, noderiv=noderiv
    n = n_elements(a)
    y = x*0.
    dyda = dblarr(n, n_elements(x))
    for l=0,n-1 do begin
        for k=0,l do dyda[l,*] += 2.^(-1.*l)*factorial(l)/(factorial(k)*factorial(l-k))^2.*(x-1.)^double(l-k)*(x+1.)^double(k)
        y += a[l]*dyda[l,*]
    endfor
    return, keyword_set(noderiv) ? y : [transpose(reform(y)), dyda]
end
