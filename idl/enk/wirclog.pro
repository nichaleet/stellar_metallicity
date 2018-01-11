pro wirclog
    f = file_search('wirc*.fits*', count=c)
    if file_test('wirclog') then begin
        expnum = '    '
        openr, 1, 'wirclog'
        skip_lun, 1, 1, /lines
        while ~eof(1) do begin
            readf, 1, expnum, format='(A4)'
            expno = fix(expnum)
            if expno ne 0 then expn = (size(expn))[1] gt 0 ? [expn, expno] : expno
        endwhile
        close, 1
        start = where(fix(strmid(f, 4, 4)) eq max(expn))
        if start ge max(fix(strmid(f, 4, 4)))-1 then begin
            print, 'No updates.'
            return
        endif else start += 1
    endif else start = 0

    openw, 1, 'wirclog', append=file_test('wirclog')
    if start eq 0 then printf, 1, 'exp# UTC      exptime coad RA          Dec         arms filter         object                                  comment'
    for i=start[0],c-1 do begin
        hdr = headfits(f[i])
        expnum = strmid(f[i], 4, 4)
        utc = strmid(sxpar(hdr, 'UTSHUT'), 11)
        coadds = sxpar(hdr, 'COADDS')
        ra = sxpar(hdr, 'RA')
        dec = sxpar(hdr, 'dec')
        airmass = sxpar(hdr, 'AIRMASS')
        exptime = sxpar(hdr, 'EXPTIME')
        object = sxpar(hdr, 'OBJECT')
        fore_filter = sxpar(hdr, 'FORE')
        aft_filter = sxpar(hdr, 'AFT')
        if strtrim(fore_filter, 2) eq 'OPEN' then filter = strtrim(aft_filter, 2)
        if strtrim(aft_filter, 2) eq 'OPEN' then filter = strtrim(fore_filter, 2)
        if strtrim(fore_filter, 2) ne 'OPEN' and strtrim(aft_filter, 2) ne 'OPEN' then filter = 'BLOCK'
        printf, 1, expnum, strmid(utc, 0, 8), exptime, coadds, ra, dec, airmass, filter, object, format='(I04,1X,A8,1X,D7.2,1X,I4,2(1X,A11),1X,D4.2,1X,A-14,1X,A-40)'
    endfor
    close, 1
end
