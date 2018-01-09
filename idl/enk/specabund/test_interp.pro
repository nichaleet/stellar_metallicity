pro test_interp
    common mooglambdacom, mooglambda
    common hashtable, ht, nmoog 
    mooglambda = read_moog_lambda7(n=nmoog)
    ht = obj_new('hashtable', length=4096, /no_duplicates, null_value=-1)

    n = 201
    alphafe = dindgen(n)*0.01-0.8
    moogspec = dblarr(20001, n)
    for i=0,n-1 do begin
        moogspec[*,i] = interp_moog7(4600, 1.5, -1.5, alphafe[i])
    endfor
    junk = min(abs(mooglambda - 8542), wl)
    splot, alphafe, moogspec[wl,*]
    w = where(round(alphafe*100) mod 10 eq 0)
    soplot, alphafe[w], moogspec[wl,w], psym=4
end
