function na_ew, spec1dfile, z, err=err, stop=stop
    spec1d = readspec(spec1dfile)
    lam_na = [8179.0, 8200.0]
    lam_cob = [8130.0, 8175.0]
    lam_cor = [8210.0, 8220.0]
    lambda = spec1d.lambda/(1d + z)
    wna = where(lambda ge lam_na[0] and lambda le lam_na[1], cna)
    wcob = where(lambda ge lam_cob[0] and lambda le lam_cob[1], ccob)
    wcor = where(lambda ge lam_cor[0] and lambda le lam_cor[1], ccor)

    if cna eq 0 or ccob eq 0 or ccor eq 0 then begin
        err = -999d
        return, -999d
    endif

    lamb = total(lambda[wcob]*spec1d.ivar[wcob])/total(spec1d.ivar[wcob])
    lamr = total(lambda[wcor]*spec1d.ivar[wcor])/total(spec1d.ivar[wcor])

    nmc = 1000
    ew_na = dblarr(nmc)
    seed = long(systime(1))
    for i=0,nmc-1 do begin
        specb = spec1d.spec[wcob] + (spec1d.ivar[wcob])^(-0.5)*randomn(seed, ccob)
        ;specb = spec1d.spec[wcob]
        contb = total(specb*spec1d.ivar[wcob])/total(spec1d.ivar[wcob])
        specr = spec1d.spec[wcor] + (spec1d.ivar[wcor])^(-0.5)*randomn(seed, ccor)
        ;specr = spec1d.spec[wcor]
        contr = total(specr*spec1d.ivar[wcor])/total(spec1d.ivar[wcor])
        slope = (contr-contb)/(lamr-lamb)
        cont = slope*(lambda[wna]-lamb)+contb
        spec = spec1d.spec[wna] + (spec1d.ivar[wna])^(-0.5)*randomn(seed, cna)
        ;spec = spec1d.spec[wna]
        wnas = [wna, wna[cna-1]+1]
        dlam = shift(lambda[wnas],-1) - lambda[wnas]
        dlam = dlam[0:cna-1]
        ew_na[i] = total((cont-spec)*spec1d.ivar[wna]*dlam/cont)*cna/total(spec1d.ivar[wna])
        ;ew_na[i] = total((cont-spec)*dlam/cont)
        ;if keyword_set(stop) then begin
        ;    plot, lambda[wna], spec1d.spec[wna]
        ;    oplot, lambda[wna], cont
        ;    stop
        ;endif
    endfor
    err = stddev(ew_na)
    ew_na = mean(ew_na)
    return, ew_na
end
