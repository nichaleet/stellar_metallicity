pro ascii_hires, inname
    if (size(inname))[1] eq 0 then message, 'You must specify an input filename.'

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;;;;;;  READ SPECTRUM  ;;;;;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    openw, lun, inname+'.dat', /get_lun

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
            printf, lun, '#chip '+strtrim(i+1, 2)+', order '+strtrim(j+1, 2)
            for k=0,npix-1 do begin
                printf, lun, lambda[k,j], flux[k,j]
            endfor
        endfor
    endfor

    close, lun
    free_lun, lun
end
