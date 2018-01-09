pro merge_hires, star, mc=mc, daospec=daospec;, minvr=minvr, maxvr=maxvr
    if ~keyword_set(mc) then mc = 0
    ;if ~keyword_set(minvr) then minvr = -350
    ;if ~keyword_set(maxvr) then maxvr = 350

    restore, getenv('CALTECH')+'hires/'+star+'_vr.sav'
    minvr = vr - 0.5
    maxvr = vr + 0.5

    for i=1,3 do begin
        starfile = star+'_'+string(i, format='(I1)')+'.fits'
        errfile = star+'_'+string(i, format='(I1)')+'.fits'
        flux = mrdfits(starfile, 0, hdrtemp, /silent)
        if i eq 1 then begin
            hdrs = hdrtemp
            err = mrdfits(errfile, 0, hdre, /silent)
        endif
        lambda = i eq 1 ? hires_lambda(hdrtemp) : [[lambda], [hires_lambda(hdrtemp)]]
    endfor
    npix = (size(lambda))[1]
    norders = (size(lambda))[2]
    minlambda = min(lambda)
    maxlambda = max(lambda)
    for i=0,norders-1 do begin
        dlambdatemp = -1.0*ts_diff(lambda[*,i], 1)
        rtemp = lambda[*,i] / dlambdatemp
        rtemp = rtemp[0:npix-2]
        dlambdatemp = dlambdatemp[0:npix-2]
        dlambda = i eq 0 ? dlambdatemp : [dlambda, dlambdatemp]
        r = i eq 0 ? rtemp : [r, rtemp]
    endfor
    dl = min(dlambda)
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

    sxaddpar, hdrs, 'CRPIX1', 1, 'Modified by merge_hires.pro'
    sxaddpar, hdrs, 'CRVAL1', min(lambda), 'Modified by merge_hires.pro'
    sxaddpar, hdrs, 'CDELT1', dl, 'Modified by merge_hires.pro'
    sxaddpar, hdrs, 'CTYPE1', 'LINEAR  ', 'Modified by merge_hires.pro'
    sxaddpar, hdrs, 'WCSDIM', 1, 'Modified by merge_hires.pro'
    sxaddpar, hdrs, 'CD1_1', dl, 'Modified by merge_hires.pro'
    sxaddpar, hdrs, 'LTM1_1', 1., 'Modified by merge_hires.pro'
    sxaddpar, hdrs, 'WAT0_001', 'system=equispec', 'Modified by merge_hires.pro'
    sxaddpar, hdrs, 'WAT1_001', 'wtype=linear label=Wavelength units=angstroms', 'Modified by merge_hires.pro'
    sxdelpar, hdrs, 'NAXIS2'

    sxaddpar, hdre, 'CRPIX1', 1, 'Modified by merge_hires.pro'
    sxaddpar, hdre, 'CRVAL1', min(lambda), 'Modified by merge_hires.pro'
    sxaddpar, hdre, 'CDELT1', dl, 'Modified by merge_hires.pro'
    sxaddpar, hdre, 'CTYPE1', 'LINEAR  ', 'Modified by merge_hires.pro'
    sxaddpar, hdre, 'WCSDIM', 1, 'Modified by merge_hires.pro'
    sxaddpar, hdre, 'CD1_1', dl, 'Modified by merge_hires.pro'
    sxaddpar, hdre, 'LTM1_1', 1., 'Modified by merge_hires.pro'
    sxaddpar, hdre, 'WAT0_001', 'system=equispec', 'Modified by merge_hires.pro'
    sxaddpar, hdre, 'WAT1_001', 'wtype=linear label=Wavelength units=angstroms', 'Modified by merge_hires.pro'
    sxdelpar, hdre, 'NAXIS2'

    mwrfits, float(spec), (mc gt 0 ? string(mc, format='(I05)')+'/' : '')+star+'_merged.fits', hdrs, /create
    mwrfits, float(resid), (mc gt 0 ? string(mc, format='(I05)')+'/' : '')+star+'_resid_merged.fits', hdrs, /create
    mwrfits, float(err), (mc gt 0 ? string(mc, format='(I05)')+'/' : '')+star+'_err_merged.fits', hdre, /create

    if keyword_set(daospec) then begin
        case strtrim(sxpar(hdrs, 'DECKNAME'), 2) of
            'C1': R = 50000d
            'C2': R = 50000d
            'C5': R = 37500d
            'E4': R = 86600d
            else: R = 37500d
        endcase
        midlambda = lambda[floor(double(n_elements(lambda))/2.)]
        ;fwhm = (midlambda / R) / dl * 2d * sqrt(2d * alog(2d))
        fwhm = 6.5

        openw, lun, (mc gt 0 ? string(mc, format='(I05)')+'/' : '')+'daospec_merged.opt', /get_lun
        printf, lun, 'or=20'
        printf, lun, 'fw='+string(fwhm, format='(D4.1)')
        printf, lun, 'sh='+string(ceil(min(lambda)), format='(I4)')
        printf, lun, 'lo='+string(floor(max(lambda)), format='(I4)')
        printf, lun, 'le='+string(floor(midlambda-5.0), format='(I4)')
        printf, lun, 'ri='+string(floor(midlambda+5.0), format='(I4)')
        printf, lun, 're = 5'
        printf, lun, 'mi = '+strtrim(string(minvr, format='(D10.3)'), 2)
        printf, lun, 'ma = '+strtrim(string(maxvr, format='(D10.3)'), 2)
        printf, lun, 've=3'
        printf, lun, 'fi=0.'
        printf, lun, 'cr=1'
        printf, lun, 'wa=1'
        printf, lun, 'sm=5'
        printf, lun, 'sc=1.'
        printf, lun, 'ba=-100000'
        close, lun
        free_lun, lun
    endif

    splot, lambda / (1d + vr / 3d5), spec, xrange=[5900, 5935], xstyle=1, yrange=[0.0, 1.3], ystyle=1
    soplot, [0, 10000], [1.0, 1.0], color=fsc_color('red')
end
