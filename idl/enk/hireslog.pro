pro hireslog
    f = file_search('hires*.fits*', count=c)
    if file_test('hireslog') then begin
        expnum = '    '
        openr, 1, 'hireslog'
        skip_lun, 1, 1, /lines
        while ~eof(1) do begin
            readf, 1, expnum, format='(A4)'
            expno = fix(expnum)
            if expno ne 0 then expn = (size(expn))[1] gt 0 ? [expn, expno] : expno
        endwhile
        close, 1
        start = where(fix(strmid(f, 5, 4)) eq max(expn))
        if start ge max(fix(strmid(f, 5, 4)))-1 then begin
            print, 'No updates.'
            return
        endif else start += 1
    endif else start = 0

    openw, 1, 'hireslog', append=file_test('hireslog')
    if start eq 0 then printf, 1, 'exp# UTC      etim xdangl ecangl filtr lamp    lfl dc arms object'
    for i=start[0],c-1 do begin
        hdr = headfits(f[i])
        expnum = strmid(f[i], 5, 4)
        utc = sxpar(hdr, 'UTC')  ;11
        xdangle = sxpar(hdr, 'XDANGL')
        airmass = sxpar(hdr, 'AIRMASS')
        exptime = sxpar(hdr, 'EXPTIME')
        object = sxpar(hdr, 'OBJECT')
        deckname = sxpar(hdr, 'DECKNAME')
        ecangle = sxpar(hdr, 'ECHANGL')
        filter = sxpar(hdr, 'FIL1NAME')
        lamp = sxpar(hdr, 'LAMPNAME')
        lampfilter = sxpar(hdr, 'LFILNAME')
        ; = sxpar(hdr, '')
        ; = sxpar(hdr, '')
        printf, 1, expnum, strmid(utc, 0, 8), exptime, xdangle, ecangle, strmid(filter, 0, 5), strmid(lamp, 0, 7), strmid(lampfilter, 0, 3), strmid(deckname, 0, 2), airmass, object, format='(I04,1X,A8,1X,I4,1X,2(D+6.3,1X),A5,1X,A7,1X,A3,1X,A2,1X,D4.2,1X,A-20)'
    endfor
    close, 1
end
