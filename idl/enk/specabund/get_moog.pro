function get_moog, xin, a
    common moog, datalambda, dlam, contmask, contdivivar, smoothfactor, phot

    teff = a[0]
    logg = a[1]
    feh = a[2]
    alpha = a[3]
    resscale = a[4]

    if phot eq 1 then x = xin[2:n_elements(xin)-1] else x = xin

    moogspec = dblarr(n_elements(x)) - 99999.
    if teff ge 7000 and logg lt 0.5 then return, phot eq 1 ? [teff, logg, moogspec] : moogspec

    moogspec = interp_moog(teff, logg, feh, alpha, mooglambda=mooglambda)
    if 0 then begin
        w1 = where(mooglambda gt min(x)-10 and mooglambda lt max(x)+10)
        w2 = where(datalambda gt min(x)-10 and datalambda lt max(x)+10)

        dl = datalambda[w2]
        contmaski = contmask[w2]
        contdivivari = contdivivar[w2]

        mooglambda = mooglambda[w1]
        moogspec = smooth_gauss(mooglambda, moogspec[w1], dl, dlam[w2]*resscale)
    endif else begin
        dl = datalambda
        contmaski = contmask
        contdivivari = contdivivar
        w2 = lindgen(n_elements(dl))

        mooglambda = mooglambda
        moogspec = smooth_gauss(mooglambda, moogspec, dl, dlam*resscale)
    endelse

    n = n_elements(datalambda)
    ;moogspecc = moogspec ;+ (contdivivar)^(-0.5)*randomn(seed, n, /double, /normal)
    nhalf = round(double(n) / 2.)
    ;wwhole = lindgen(nhalf)
    for ccd=1,2 do begin
        case ccd of
            1: begin
                ;ww = where(w2 ge 3 and w2 lt nhalf-3)
                wwhole = where(w2 lt nhalf)
            end
            2: begin
                ;ww = where(w2 ge nhalf+3 and w2 lt n-3)
                wwhole = where(w2 ge nhalf)
            end
        endcase
;        wcont = ww[where(moogspec[ww] gt 0 and dl[ww] gt min(mooglambda) and dl[ww] lt max(mooglambda) and contmaski[ww] eq 1 and finite(contdivivari[ww]) and contdivivari[ww] gt 0)]
;        for i=1,5 do begin
;            wcontold = wcont
;            cont = smooth_gauss(dl[wcont], moogspec[wcont], dl[wwhole], 30.0, ivar1=contdivivari[wcont])
;            wcont = ww[where(moogspec[ww] gt 0 and dl[ww] gt min(mooglambda) and dl[ww] lt max(mooglambda) and abs(moogspec[ww]-cont[ww-wwhole[0]]) lt (contdivivari[ww])^(-0.5) and finite(contdivivari[ww]) and contdivivari[ww] gt 0 and contmaski[ww] eq 1)]
;            if array_equal(wcont, wcontold) then break
;        endfor

        won = where(contmaski[wwhole] eq 1 and finite(contdivivar[wwhole]) and contdivivar[wwhole] gt 0.0, complement=woff, con) + wwhole[0]
        bkpt = slatec_splinefit(dl[won], moogspec[won], coeff, invvar=contdivivar[won], bkspace=150, upper=10, lower=10, /silent)
        cont = slatec_bvalu(dl[wwhole], bkpt, coeff)
        case ccd of
            1: contb = cont
            2: contr = cont
        endcase
    endfor
    moogcont = [contb, contr]
    moogspec /= moogcont

    if smoothfactor gt 0 then moogspec = smooth_gauss(dl, moogspec, x, smoothfactor, ivar1=contdivivar)
    match, dl, x, w, wjunk
    moogspec = moogspec[w]
    return, phot eq 1 ? [teff, logg, moogspec] : moogspec
end
