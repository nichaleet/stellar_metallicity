pro dbsplog
    bf = file_search('blue????.fits*', count=cb)
    if file_test('dbsp.blue.log') then begin
        expnum = '    '
        openr, 1, 'dbsp.blue.log'
        skip_lun, 1, 4, /lines
        while ~eof(1) do begin
            readf, 1, expnum, format='(A4)'
            expno = fix(expnum)
            if expno ne 0 then expn = (size(expn))[1] gt 0 ? [expn, expno] : expno
        endwhile
        close, 1
        start = where(fix(strmid(bf, 4, 4)) eq max(expn))
        if start ge max(fix(strmid(bf, 4, 4)))-1 then begin
            print, 'No blue updates.'
            goto, red
        endif else start += 1
    endif else start = 0

    openw, 1, 'dbsp.blue.log', append=file_test('dbsp.blue.log')
    for i=start[0],cb-1 do begin
        hdr = headfits(bf[i])
        expnum = strmid(bf[i], 4, 4)
        grating = sxpar(hdr, 'GRATING')
        angle = sxpar(hdr, 'ANGLE')
        utc = sxpar(hdr, 'UT')
        ra = sxpar(hdr, 'RA')
        dec = sxpar(hdr, 'DEC')
        casspa = sxpar(hdr, 'CASSPA')
        aperture = sxpar(hdr, 'APERTURE')
        airmass = sxpar(hdr, 'AIRMASS')
        exptime = sxpar(hdr, 'EXPTIME')
        object = sxpar(hdr, 'OBJECT')
        lamps = sxpar(hdr, 'LAMPS')   ;D-FeAr-Hg-Ar-Ne-He-InCan
        if i eq 0 then begin
            printf, 1, 'DBSP blue'
            printf, 1, 'GRATING: ', grating, 'ANGLE: ', angle, format='(A9,A9,5X,A7,A13)'
            printf, 1, ' '
            printf, 1, 'exp# UTC      exptime RA          Dec         casspa arms aper lamps   object                                  comment'
        endif
        printf, 1, expnum, strmid(utc, 0, 8), exptime, ra, dec, casspa, airmass, aperture, lamps, object, format='(I04,1X,A8,2X,D6.1,2(1X,A11),1X,D6.1,1X,D4.2,1X,A4,1X,A7,1X,A-40)'
    endfor
    close, 1

    delvarx, expn
    red:
    rf = file_search('red????.fits*', count=cr)
    if file_test('dbsp.red.log') then begin
        expnum = '    '
        openr, 1, 'dbsp.red.log'
        skip_lun, 1, 4, /lines
        while ~eof(1) do begin
            readf, 1, expnum, format='(A4)'
            expno = fix(expnum)
            if expno ne 0 then expn = (size(expn))[1] gt 0 ? [expn, expno] : expno
        endwhile
        close, 1
        start = where(fix(strmid(rf, 3, 4)) eq max(expn))
        if start ge max(fix(strmid(rf, 3, 4)))-1 then begin
            print, 'No red updates.'
            return
        endif else start += 1
    endif else start = 0

    openw, 1, 'dbsp.red.log', append=file_test('dbsp.red.log')
    for i=start[0],cr-1 do begin
        hdr = headfits(rf[i])
        expnum = strmid(rf[i], 3, 4)
        grating = sxpar(hdr, 'GRATING')
        angle = sxpar(hdr, 'ANGLE')
        utc = sxpar(hdr, 'UT')
        ra = sxpar(hdr, 'RA')
        dec = sxpar(hdr, 'DEC')
        casspa = sxpar(hdr, 'CASSPA')
        aperture = sxpar(hdr, 'APERTURE')
        airmass = sxpar(hdr, 'AIRMASS')
        exptime = sxpar(hdr, 'EXPTIME')
        object = sxpar(hdr, 'OBJECT')
        lamps = sxpar(hdr, 'LAMPS')   ;D-FeAr-Hg-Ar-Ne-He-InCan
        if i eq 0 then begin
            printf, 1, 'DBSP red'
            printf, 1, 'GRATING: ', grating, 'ANGLE: ', angle, format='(A9,A9,5X,A7,A13)'
            printf, 1, ' '
            printf, 1, 'exp# UTC      exptime RA          Dec         casspa arms aper lamps   object                                  comment'
        endif
        printf, 1, expnum, strmid(utc, 0, 8), exptime, ra, dec, casspa, airmass, aperture, lamps, object, format='(I04,1X,A8,2X,D6.1,2(1X,A11),1X,D6.1,1X,D4.2,1X,A4,1X,A7,1X,A-40)'
    endfor
    close, 1
end
