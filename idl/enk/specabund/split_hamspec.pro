pro split_hamspec, inname, outname=outname, daospec=daospec, mcnum=mcnum, seed=seed, minvr=minvr, maxvr=maxvr
    if (size(inname))[1] eq 0 then message, 'You must specify an input filename.'
    if ~keyword_set(outname) then outname = inname
    if ~keyword_set(minvr) then minvr = -350
    if ~keyword_set(maxvr) then maxvr = 350
    if keyword_set(mcnum) and ~keyword_set(seed) then message, 'If you specify a Monte Carlo trial number, then you need to specify the random number seed.'
    if keyword_set(mcnum) then begin
        outname = string(mcnum, format='(I05)')+'/'+outname
        file_delete, string(mcnum, format='(I05)'), /recursive, /quiet
        file_mkdir, string(mcnum, format='(I05)')
    endif

    case inname of
        'HIP57939a': vr = -14.0
        else:
    endcase
    if keyword_set(vr) then begin
        minvr = vr - 0.5
        maxvr = vr + 0.5
    endif

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;;;;;;  READ SPECTRUM  ;;;;;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    orderi = 0
    data = mrdfits(inname+'.fits.gz', 1, /silent)
    lambda = data.wave
    flux = data.fx
    err = sqrt(data.var)
    norders = (size(flux))[2]

    for j=0,norders-1 do begin
        orderi++
        lambdaj = lambda[*,j]
        w = where(lambdaj gt 0, npix)
        if npix eq 0 then continue
        lambdaj = lambdaj[w]
        fluxj = flux[w,j]
        errj = err[w,j]
        mkhdr, hdrf, fluxj
        dlambda = (max(lambdaj)-min(lambdaj)) / double(npix)
        lambdalin = dindgen(npix)*dlambda + min(lambdaj)
        x_specrebin, lambdaj, fluxj, lambdalin, fluxlin, /silent
        x_specrebin, lambdaj, errj, lambdalin, errlin, /silent
        if keyword_set(mcnum) then begin
            fluxlin += errlin*randomn(seed, n_elements(fluxlin))
        endif
        w = where(~finite(fluxlin) or ~finite(errlin) or errlin le 0.0 or (lambda gt 6866 and lambda lt 6969) or (lambda gt 7591 and lambda lt 7681))

        sxaddpar, hdrf, 'CRPIX1', 1, 'Modified by split_hamspec.pro'
        sxaddpar, hdrf, 'CRVAL1', min(lambdaj), 'Modified by split_hamspec.pro'
        sxaddpar, hdrf, 'CDELT1', dlambda, 'Modified by split_hamspec.pro'
        sxaddpar, hdrf, 'CTYPE1', 'LINEAR  ', 'Modified by split_hamspec.pro'
        sxaddpar, hdrf, 'WCSDIM', 1, 'Modified by split_hamspec.pro'
        sxaddpar, hdrf, 'CD1_1', dlambda, 'Modified by split_hamspec.pro'
        sxaddpar, hdrf, 'LTM1_1', 1., 'Modified by split_hamspec.pro'
        sxaddpar, hdrf, 'WAT0_001', 'system=equispec', 'Modified by split_hamspec.pro'
        sxaddpar, hdrf, 'WAT1_001', 'wtype=linear label=Wavelength units=angstroms', 'Modified by split_hamspec.pro'
        mwrfits, fluxlin, outname+'_'+string(orderi, format='(I03)')+'.fits', hdrf, /create
        if ~keyword_set(mcnum) then mwrfits, errlin, outname+'_'+string(orderi, format='(I03)')+'e.fits', hdrf, /create

        if keyword_set(daospec) then begin
            R = 60000d
            midlambda = lambdalin[floor(double(npix)/2.)]
            ;fwhm = (midlambda / R) / dlambda * 2d * sqrt(2d * alog(2d))
            fwhm = 3.5

            openw, lun, (keyword_set(mcnum) ? string(mcnum, format='(I05)')+'/' : '')+'daospec_'+string(orderi, format='(I03)')+'.opt', /get_lun
            printf, lun, 'or=20'
            printf, lun, 'fw='+string(fwhm, format='(D4.1)')
            printf, lun, 'sh='+string(ceil(min(lambdalin)), format='(I4)')
            printf, lun, 'lo='+string(floor(max(lambdalin)), format='(I4)')
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
    endfor
end
