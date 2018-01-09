pro moog_add_err, moogify, aherr=aherr, nocovar=nocovar
    tagnames = tag_names(moogify)

    if 0 and contains(tagnames, 'COVAR') and ~keyword_set(nocovar) then begin
        n = n_elements(moogify)
        feherr = dblarr(n)
        tefferr = dblarr(n)
        for i=0,n-1 do begin
            covar = moogify[i].covar
            covar = [[covar[0,0], covar[0,2]], [covar[2,0], [covar[2,2]]]]
            alpha = invert(covar, /double)
            teffextremum = abs(alpha[1,0]*sqrt(2.30/(alpha[0,0]*(alpha[0,0]*alpha[1,1]-alpha[1,0]^2.))))
            feherr[i] = (-1.*alpha[1,0]*teffextremum + sqrt((alpha[1,0]*teffextremum)^2.-alpha[1,1]*(alpha[0,0]*teffextremum^2.-2.30)))/alpha[1,1]
            fehextremum = abs(alpha[1,0]*sqrt(2.30/(alpha[1,1]*(alpha[0,0]*alpha[1,1]-alpha[1,0]^2.))))
            tefferr[i] = (-1.*alpha[1,0]*fehextremum + sqrt((alpha[1,0]*fehextremum)^2.-alpha[0,0]*(alpha[1,1]*fehextremum^2.-2.30)))/alpha[0,0]
            ;resetplot
            ;dfeh = (dindgen(1001)-500.)/500.
            ;dteff1 = (-1.*alpha[1,0]*dfeh - sqrt((alpha[1,0]*dfeh)^2.-alpha[0,0]*(alpha[1,1]*dfeh^2.-2.30)))/alpha[0,0]
            ;dteff2 = (-1.*alpha[1,0]*dfeh + sqrt((alpha[1,0]*dfeh)^2.-alpha[0,0]*(alpha[1,1]*dfeh^2.-2.30)))/alpha[0,0]
            ;plot, dfeh, dteff1, xrange=[-1.*feherr[i], feherr[i]], xstyle=1
            ;oplot, dfeh, dteff2
            ;stop
        endfor
    
        moogify.tefferr = tefferr
        moogify.feherr = feherr
    endif

    ;fehsyserr = 0.11334408
    ;alphasyserr = 0.081761405
    ;mgfesyserr = 0.095349035
    ;sifesyserr = 0.10635760
    ;cafesyserr = 0.11795958
    ;tifesyserr = 0.083345975
    fehsyserr = 0.10630998         ;updated 11/08/12 after fixing a bug in interp_moog8.pro
    alphasyserr = 0.086330372
    mgfesyserr = 0.064957422
    sifesyserr = 0.11339209
    cafesyserr = 0.11104078
    tifesyserr = 0.090354729
    ;fehsyserr = alog(exp(0.1000000) + exp(-0.8945806 + 1.9643203*moogify.feh))
    w = where(moogify.feherr gt 0.0, c)
    if c gt 0 then moogify[w].feherr = sqrt(moogify[w].feherr^2. + fehsyserr^2.)
    w = where(moogify.alphafeerr gt 0.0, c)
    if c gt 0 then moogify[w].alphafeerr = sqrt(moogify[w].alphafeerr^2. + alphasyserr^2.)
    if contains(tagnames, 'MGFE') then begin
        w = where(moogify.mgfeerr gt 0.0, c)
        if c gt 0 then moogify[w].mgfeerr = sqrt(moogify[w].mgfeerr^2. + mgfesyserr^2.)
    endif
    if contains(tagnames, 'SIFE') then begin
        w = where(moogify.sifeerr gt 0.0, c)
        if c gt 0 then moogify[w].sifeerr = sqrt(moogify[w].sifeerr^2. + sifesyserr^2.)
    endif
    if contains(tagnames, 'CAFE') then begin
        w = where(moogify.cafeerr gt 0.0, c)
        if c gt 0 then moogify[w].cafeerr = sqrt(moogify[w].cafeerr^2. + cafesyserr^2.)
    endif
    if contains(tagnames, 'TIFE') then begin
        w = where(moogify.tifeerr gt 0.0, c)
        if c gt 0 then moogify[w].tifeerr = sqrt(moogify[w].tifeerr^2. + tifesyserr^2.)
    endif

    tags = tag_names(moogify)
    if contains(tags, 'COVAR') then begin
        aherr = sqrt(moogify.covar[2,2] + moogify.covar[3,3] + moogify.covar[2,3] + moogify.covar[3,2])*sqrt(moogify.chisq)
        aherr = sqrt(aherr^2. + fehsyserr^2.)
    endif
end

