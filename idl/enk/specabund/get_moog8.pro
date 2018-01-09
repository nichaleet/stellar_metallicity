function get_moog8, xin, a, addtofeh=addtofeh
    common moog, datalambda, dlam, moogmask, contmask, contdiv, contdivivar, smoothfactor, phot, moogcontrevise, pixmask, hectospec
    common mooglambdacom, mooglambda

    teff = a[0]
    logg = a[1]
    feh = a[2]
    alpha = a[3]
    mgsicati = a[4]
    resscale = a[5]

    if keyword_set(addtofeh) then feh += addtofeh

    if phot eq 1 then x = xin[2:n_elements(xin)-1] else x = xin

    moogspec = dblarr(n_elements(x)) - 99999.
    if teff ge 6000 and logg lt 0.5 then return, phot eq 1 ? [teff, logg, moogspec] : moogspec
    if logg ge 3.6 and feh lt -4.0 then return, phot eq 1 ? [teff, logg, moogspec] : moogspec

    moogspec = interp_moog8(teff, logg, feh, alpha, mgsicati)
    ;moogspec = interp_moog8(teff, logg, feh, alpha, mgsicati, mooglambda=mooglambda, /fullres)
    if 0 then begin
        w1 = where(mooglambda gt min(x)-10 and mooglambda lt max(x)+10)
        w2 = where(datalambda gt min(x)-10 and datalambda lt max(x)+10)

        dl = datalambda[w2]
        contmaski = contmask[w2]
        contdivi = contdiv[w2]
        contdivivari = contdivivar[w2]

        mooglambda = mooglambda[w1]
        moogspec = smooth_gauss_wrapper(mooglambda, moogspec[w1], dl, dlam[w2]*resscale)
    endif else begin
        dl = datalambda
        contmaski = contmask
        contdivi = contdiv
        contdivivari = contdivivar
        w2 = lindgen(n_elements(dl))

        ;mooglambda = mooglambda
        moogspec = smooth_gauss_wrapper(mooglambda, moogspec, dl, dlam*resscale)
    endelse

    n = n_elements(dl)
    if 0 and moogcontrevise eq 1 then begin
        ;moogspecc = moogspec ;+ (contdivivar)^(-0.5)*randomn(seed, n, /double, /normal)
        nhalf = hectospec eq 1 ? n : round(double(n) / 2.)
        ;wwhole = lindgen(nhalf)
        for ccd=1,(hectospec eq 1 ? 1 : 2) do begin
            case ccd of
                1: begin
                    ;ww = where(w2 ge 3 and w2 lt nhalf-3)
                    wwhole = where(w2 ge 5 and w2 lt nhalf-5)
                end
                2: begin
                    ;ww = where(w2 ge nhalf+3 and w2 lt n-3)
                    wwhole = where(w2 ge nhalf+5 and w2 lt n)
                end
            endcase
;        wcont = ww[where(moogspec[ww] gt 0 and dl[ww] gt min(mooglambda) and dl[ww] lt max(mooglambda) and contmaski[ww] eq 1 and finite(contdivivari[ww]) and contdivivari[ww] gt 0)]
;        for i=1,5 do begin
;            wcontold = wcont
;            cont = smooth_gauss_wrapper(dl[wcont], moogspec[wcont], dl[wwhole], 30.0, ivar1=contdivivari[wcont])
;            wcont = ww[where(moogspec[ww] gt 0 and dl[ww] gt min(mooglambda) and dl[ww] lt max(mooglambda) and abs(moogspec[ww]-cont[ww-wwhole[0]]) lt (contdivivari[ww])^(-0.5) and finite(contdivivari[ww]) and contdivivari[ww] gt 0 and contmaski[ww] eq 1)]
;            if array_equal(wcont, wcontold) then break
;        endfor
    
            won = where(contmaski[wwhole] eq 1 and finite(contdivivari[wwhole]) and contdivivari[wwhole] gt 0.0, complement=woff, con) + wwhole[0]
            bkpt = slatec_splinefit(dl[won], moogspec[won], coeff, invvar=contdivivari[won], bkspace=150, upper=5, lower=5, /silent)
            cont = slatec_bvalu(dl[wwhole], bkpt, coeff)
            case ccd of
                1: contb = cont
                2: contr = cont
            endcase
        endfor
        moogcont = hectospec eq 1 ? contb : [contb, contr]
        moogspec /= moogcont
    endif

    if 0 then begin
        w2 = lindgen(n)
        nhalf = hectospec eq 1 ? n : round(double(n) / 2.)
        for ccd=1,(hectospec eq 1 ? 1 : 2) do begin
            case ccd of
                1: wwhole = where(w2 ge 5 and w2 lt nhalf-5)
                2: wwhole = where(w2 ge nhalf+5 and w2 lt n)
            endcase
            won = where(moogspec[wwhole] gt 0 and pixmask[wwhole] eq 1 and finite(contdivivar[wwhole]) and contdivivar[wwhole] gt 0.0, complement=woff, con) + wwhole[0]
            if con eq 0 then cont = dblarr(n_elements(wwhole))+1.0 else begin 
                bkpt = slatec_splinefit(dl[won], contdiv[won]/moogspec[won], coeff, invvar=contdivivar[won], bkspace=25, upper=5, lower=5, /silent)
                cont = slatec_bvalu(dl[wwhole], bkpt, coeff)
            endelse
            case ccd of
                1: contb = cont
                2: contr = cont
            endcase
        endfor
        speccont = hectospec eq 1 ? contb : [contb, contr]
        moogspec *= speccont
    endif

    ;for ccd=1,2 do begin
    ;    case ccd of
    ;        1: begin
    ;            ;ww = where(w2 ge 3 and w2 lt nhalf-3)
    ;            wwhole = where(w2 lt nhalf)
    ;        end
    ;        2: begin
    ;            ;ww = where(w2 ge nhalf+3 and w2 lt n-3)
    ;            wwhole = where(w2 ge nhalf)
    ;        end
    ;    endcase
    ;    won = where(moogmask[wwhole] eq 1 and finite(contdivivar[wwhole]) and contdivivar[wwhole] gt 0.0, complement=woff, con) + wwhole[0]
    ;    bkpt = slatec_splinefit(dl[won], contdivi[won]/moogspec[won], coeff, invvar=contdivivar[won], bkspace=50, upper=5, lower=5, /silent)
    ;    cont = slatec_bvalu(dl[wwhole], bkpt, coeff)
    ;    case ccd of
    ;        1: contb = cont
    ;        2: contr = cont
    ;    endcase
    ;endfor
    ;speccont = [contb, contr]
    ;moogspec *= speccont


    if mean(smoothfactor) gt 0 then moogspec = smooth_gauss_wrapper(dl, moogspec, x, smoothfactor, ivar1=contdivivar) else begin
        match, dl, x, wmask, wjunk
        moogspec = moogspec[wmask]
    endelse
    return, phot eq 1 ? [teff, logg, moogspec] : moogspec
end
