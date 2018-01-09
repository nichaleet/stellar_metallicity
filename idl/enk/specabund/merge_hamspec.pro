pro merge_hamspec, star, mc=mc
    if ~keyword_set(mc) then mc = 0

    data = mrdfits(star+'.fits.gz', 1, /silent)
    lambda = data.wave
    minlambda = min(lambda[where(lambda gt 0)])
    maxlambda = max(lambda)
    dlambda = -1.0*ts_diff(lambda[where(lambda gt 0)], 1)
    dl = min(dlambda[where(dlambda gt 0)])
    n = ceil((maxlambda-minlambda) / dl)
    lambda = dindgen(n)*dl + minlambda
    spec = dblarr(n)
    resid = dblarr(n)
    ivar = dblarr(n)

    contfiles = file_search((mc gt 0 ? string(mc, format='(I05)')+'/' : '')+star+'*C.fits', count=c)
    for i=0,c-1 do begin
        ;print, 'order '+string(i, format='(I-2)')
        fluxfile = strmid(contfiles[i], 0, strlen(contfiles[i])-6)+'.fits'
        residfile = strmid(contfiles[i], 0, strlen(contfiles[i])-6)+'R.fits'
        errfile = strmid(file_basename(contfiles[i]), 0, strlen(file_basename(contfiles[i]))-6)+'e.fits'

        fluxi = mrdfits(fluxfile, 0, hdr, /silent)
        residi = mrdfits(residfile, 0, hdr, /silent)
        erri = mrdfits(errfile, 0, /silent)
        conti = mrdfits(contfiles[i], 0, /silent)
        fluxi /= conti
        erri /= conti
        npix = sxpar(hdr, 'NAXIS1')
        dlambda = sxpar(hdr, 'CDELT1')
        lambda0 = sxpar(hdr, 'CRVAL1')
        lambdai = dindgen(npix)*dlambda+lambda0
        w = where(finite(lambdai) and finite(fluxi) and finite(erri) and erri gt 0.0, c)
        if c lt 100 then begin
            message, 'Something is wrong with order '+strtrim(string(i+1, format='(I3)'), 2)+'.', /info
            continue
        endif
        lambdai = lambdai[w]
        fluxi = fluxi[w]
        erri = erri[w]

        w = where(lambda ge min(lambdai)-5.0*dl and lambda le max(lambdai)+5.0*dl)
        x_specrebin, lambdai, fluxi, lambda[w], speci, var=(mc gt 0 ? 2.0 : 1.0)*erri^2., nwvar=nwvar, /silent
        x_specrebin, lambdai, dblarr(n_elements(lambdai))+1.0, lambda[w], norm, /silent
        ww = where(nwvar gt 0.0 and finite(speci) and finite(nwvar))
        speci[ww] /= norm[ww]
        nwvar[ww] /= (norm[ww])^2.
        spec[w[ww]] += speci[ww] / nwvar[ww]
        resid[w[ww]] += residi[ww] / nwvar[ww]
        ivar[w[ww]] += 1.0 / nwvar[ww]
    endfor
    w = where(ivar gt 0.0 and finite(spec) and finite(ivar))
    spec[w] /= ivar[w]
    resid[w] /= ivar[w]
    err = ivar
    err[w] = (ivar[w])^(-0.5)
    w = where(~finite(spec) or ~finite(ivar) or ivar le 0.0 or err gt 0.5 or (lambda gt 6410 and lambda lt 6425) or (lambda gt 7593 and lambda lt 7681))
    spec[w] = 1.0
    resid[w] = 1.0
    err[w] = 0.0

    sxaddpar, hdrs, 'CRPIX1', 1, 'Modified by merge_hamspec.pro'
    sxaddpar, hdrs, 'CRVAL1', min(lambda), 'Modified by merge_hamspec.pro'
    sxaddpar, hdrs, 'CDELT1', dl, 'Modified by merge_hamspec.pro'
    sxaddpar, hdrs, 'CTYPE1', 'LINEAR  ', 'Modified by merge_hamspec.pro'
    sxaddpar, hdrs, 'WCSDIM', 1, 'Modified by merge_hamspec.pro'
    sxaddpar, hdrs, 'CD1_1', dl, 'Modified by merge_hamspec.pro'
    sxaddpar, hdrs, 'LTM1_1', 1., 'Modified by merge_hamspec.pro'
    sxaddpar, hdrs, 'WAT0_001', 'system=equispec', 'Modified by merge_hamspec.pro'
    sxaddpar, hdrs, 'WAT1_001', 'wtype=linear label=Wavelength units=angstroms', 'Modified by merge_hamspec.pro'
    sxdelpar, hdrs, 'NAXIS2'

    sxaddpar, hdre, 'CRPIX1', 1, 'Modified by merge_hamspec.pro'
    sxaddpar, hdre, 'CRVAL1', min(lambda), 'Modified by merge_hamspec.pro'
    sxaddpar, hdre, 'CDELT1', dl, 'Modified by merge_hamspec.pro'
    sxaddpar, hdre, 'CTYPE1', 'LINEAR  ', 'Modified by merge_hamspec.pro'
    sxaddpar, hdre, 'WCSDIM', 1, 'Modified by merge_hamspec.pro'
    sxaddpar, hdre, 'CD1_1', dl, 'Modified by merge_hamspec.pro'
    sxaddpar, hdre, 'LTM1_1', 1., 'Modified by merge_hamspec.pro'
    sxaddpar, hdre, 'WAT0_001', 'system=equispec', 'Modified by merge_hamspec.pro'
    sxaddpar, hdre, 'WAT1_001', 'wtype=linear label=Wavelength units=angstroms', 'Modified by merge_hamspec.pro'
    sxdelpar, hdre, 'NAXIS2'

    mwrfits, float(spec), (mc gt 0 ? string(mc, format='(I05)')+'/' : '')+star+'_merged.fits', hdrs, /create
    mwrfits, float(resid), (mc gt 0 ? string(mc, format='(I05)')+'/' : '')+star+'_resid_merged.fits', hdrs, /create
    mwrfits, float(err), (mc gt 0 ? string(mc, format='(I05)')+'/' : '')+star+'_err_merged.fits', hdre, /create

    ;splot, lambda, spec, xrange=[5900, 5935], xstyle=1, yrange=[0.0, 1.3], ystyle=1
    ;soplot, [0, 10000], [1.0, 1.0], color=fsc_color('red')
end
