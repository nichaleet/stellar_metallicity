pro saturated_gaussian_calculate
    x = (dindgen(20001)-10000.)*0.001
    sigma = 1.0
    a = dindgen(1000001)/100d
    ew = dblarr(n_elements(a))
    for i=0L,n_elements(a)-1 do begin
        g = a[i]*exp(-(x^2d / (2d * sigma^2d)))
        y = g / (1d + g/0.95)
        ew[i] = int_tabulated(x, y, /double)
    endfor
    ;plot, a, ew, psym=-4
    save, a, ew, filename=getenv('ahome')+'idl/enk/specabund/saturated_gaussian.sav'
end
