pro split_hires, inname, outname=outname, cpsc=cpsc, daospec=daospec, mcnum=mcnum, seed=seed;, minvr=minvr, maxvr=maxvr
    if (size(inname))[1] eq 0 then message, 'You must specify an input filename.'
    if ~keyword_set(outname) then outname = inname
    ;if ~keyword_set(minvr) then minvr = -350
    ;if ~keyword_set(maxvr) then maxvr = 350
    if keyword_set(mcnum) and ~keyword_set(seed) then message, 'If you specify a Monte Carlo trial number, then you need to specify the random number seed.'
    if keyword_set(mcnum) then begin
        outname = string(mcnum, format='(I05)')+'/'+outname
        file_delete, string(mcnum, format='(I05)'), /recursive, /quiet
        file_mkdir, string(mcnum, format='(I05)')
    endif

    restore, getenv('CALTECH')+'hires/'+inname+'_vr.sav'
    minvr = vr - 0.5
    maxvr = vr + 0.5

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;;;;;;  READ SPECTRUM  ;;;;;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ntotorders = 0
    orderi = 0
    for i=0,2 do begin
        if keyword_set(cpsc) then begin
            case i of
                0: begin
                    prefix = 'b'
                    chip = 'blue'
                end
                1: begin
                    prefix = 'r'
                    chip = 'middle'
                end
                2: begin
                    prefix = 'i'
                    chip = 'red'
                end
            endcase
            flux = mrdfits(prefix+'j85.'+inname+'.fits', 0, hdrf, /silent)
            errflag = 0
        endif else begin
            flux = mrdfits(inname+'_'+string(i+1, format='(I01)')+'.fits', 0, hdrf, /silent)
            errfile = inname+'_'+string(i+1, format='(I01)')+'e.fits'
            errflag = file_test(errfile)
            if errflag then err = mrdfits(errfile, 0, /silent)
            ;errflag = 1
            ;if file_test(errfile) then err = mrdfits(errfile, 0, /silent) else err = sqrt(flux*sxpar(hdrf, 'EPERDN') + sxpar(hdrf, 'RONOISE')^2.0) / sxpar(hdrf, 'EPERDN')
        endelse
        if i eq 0 then begin
            ra = sxpar(hdrf, 'RA')
            dec = sxpar(hdrf, 'DEC')
            mjd = sxpar(hdrf, 'MJD')
        endif
        norders = (size(flux))[2]
        ;norders = sxpar(hdrf, 'MK_NEO')
        ntotorders += norders
        npix = fix(sxpar(hdrf, 'NAXIS1'))
        sxdelpar, hdrf, 'NAXIS2'
        if keyword_set(cpsc) then begin
            lambda = hires_wscale(chip=chip)
        endif else begin
            lambda = dblarr(npix, norders)
            lambdapoly = dblarr(7)
            for j=0,norders-1 do begin
                poly_string1 = sxpar(hdrf, 'WV_0_'+string(j+1, format='(I02)'))
                poly_string2 = sxpar(hdrf, 'WV_4_'+string(j+1, format='(I02)'))
                lambdapoly[0] = double(strmid(poly_string1, 0, 17))
                lambdapoly[1] = double(strmid(poly_string1, 17, 17))
                lambdapoly[2] = double(strmid(poly_string1, 34, 17))
                lambdapoly[3] = double(strmid(poly_string1, 51, 17))
                lambdapoly[4] = double(strmid(poly_string2, 0, 17))
                lambdapoly[5] = double(strmid(poly_string2, 17, 17))
                lambdapoly[6] = double(strmid(poly_string2, 34, 17))
                lambda[*,j] = poly(dindgen(npix)+1d, lambdapoly)
            endfor
        endelse

        for j=0,norders-1 do begin
            orderi++
            dlambda = (max(lambda[*,j])-min(lambda[*,j])) / double(npix)
            lambdalin = dindgen(npix)*dlambda + min(lambda[*,j])
            x_specrebin, lambda[*,j], flux[*,j], lambdalin, fluxlin, /silent
            if errflag then x_specrebin, lambda[*,j], err[*,j], lambdalin, errlin, /silent
            if keyword_set(mcnum) then begin
                fluxlin += errlin*randomn(seed, n_elements(fluxlin))
            endif
            if errflag then begin
                w = where(~finite(fluxlin) or ~finite(errlin) or errlin le 0.0 or (lambda gt 6866 and lambda lt 6969) or (lambda gt 7591 and lambda lt 7681))
            endif else begin
                w = where(~finite(fluxlin) or (lambda gt 6866 and lambda lt 6969) or (lambda gt 7591 and lambda lt 7681))
            endelse

            sxaddpar, hdrf, 'CRPIX1', 1, 'Modified by split_hires.pro'
            sxaddpar, hdrf, 'CRVAL1', min(lambda[*,j]), 'Modified by split_hires.pro'
            sxaddpar, hdrf, 'CDELT1', dlambda, 'Modified by split_hires.pro'
            sxaddpar, hdrf, 'CTYPE1', 'LINEAR  ', 'Modified by split_hires.pro'
            sxaddpar, hdrf, 'WCSDIM', 1, 'Modified by split_hires.pro'
            sxaddpar, hdrf, 'CD1_1', dlambda, 'Modified by split_hires.pro'
            sxaddpar, hdrf, 'LTM1_1', 1., 'Modified by split_hires.pro'
            sxaddpar, hdrf, 'WAT0_001', 'system=equispec', 'Modified by split_hires.pro'
            sxaddpar, hdrf, 'WAT1_001', 'wtype=linear label=Wavelength units=angstroms', 'Modified by split_hires.pro'
            mwrfits, fluxlin, outname+'_'+string(orderi, format='(I03)')+'.fits', hdrf, /create
            if ~keyword_set(mcnum) and errflag then mwrfits, errlin, outname+'_'+string(orderi, format='(I03)')+'e.fits', hdrf, /create

            if keyword_set(daospec) then begin
                case strtrim(sxpar(hdrf, 'DECKNAME'), 2) of
                    'C1': R = 50000d
                    'C2': R = 50000d
                    'C5': R = 37500d
                    'E4': R = 86600d
                    else: R = 37500d
                endcase
                midlambda = lambdalin[floor(double(npix)/2.)]
                ;fwhm = (midlambda / R) / dlambda * 2d * sqrt(2d * alog(2d))
                fwhm = 6.5

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
    endfor
end
