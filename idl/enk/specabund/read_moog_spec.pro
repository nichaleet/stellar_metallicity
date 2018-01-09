function read_moog_spec, infile, fits=fitsfile, newmoog=newmoog
    if ~file_test(infile) then message, 'File not found.'

    ;catch, error_status
    ;if error_status ne 0 then begin
    ;    message, 'Error reading '+strtrim(infile, 2), /info
    ;    close, 1
    ;    return, {elem:-1, abund:-1, lambda:-1, spec:-1}
    ;endif

    nspec = 441
    elem = strarr(95)
    abund = dblarr(nspec,95)
    openr, lun, infile, /get_lun
    temp = ' '
    teststring = ' '
    elem_in = ' '
    abund_in = 0d    
    i = 0L
    while ~eof(1) do begin
        j = 0L
        readf, 1, mh, format=(keyword_set(newmoog) ? '(53X,D6.2,4X)' : '(54X,D6.2,4X)')
        elem[j] = 'M/H'
        abund[i,j] = mh
        j += 1
        readheader:
        readf, 1, temp
        reads, temp, teststring, format='(8X,A7)'
        if teststring eq ' Ratio:' then goto, readheader
        reads, temp, teststring, elem_in, abund_in, format='(8X,A7,1X,A2,15X,D5.2)'
        while teststring eq 'element' do begin
            elem[j] = elem_in
            abund[i,j] = abund_in
            readf, 1, temp
            reads, temp, teststring, elem_in, abund_in, format='(8X,A7,1X,A2,15X,D5.2)'
            j += 1
        endwhile
        if i eq 0 then begin
            elem = elem[0:j-1]
            abund = abund[*,0:j-1]
            reads, temp, teff, logg, feh, vt, format=(keyword_set(newmoog) ? '(7X,D5.0,1X,D4.2,1X,D5.2,43X,D5.2,16X)' : '(D5.0,1X,D4.2,1X,D5.2,43X,D5.2,16X)')
            readf, 1, startlambda, endlambda, step, format='(D11.3,D11.3,D11.3)'
            startlambda = double(round(startlambda*100d)/100d)
            endlambda = double(round(endlambda*100d)/100d)
            nlambda = ceil((endlambda-startlambda)/step)
            if abs((endlambda-startlambda)/step - nlambda) lt 1d-4 then nlambda += 1
            lambda = dindgen(nlambda)*step + startlambda
            spec = fltarr(nlambda,nspec)
        endif else skip_lun, 1, 1, /lines
        for j=0L,floor(nlambda/10.0)-1 do begin
            s1 = '       '
            s2 = '       '
            s3 = '       '
            s4 = '       '
            s5 = '       '
            s6 = '       '
            s7 = '       '
            s8 = '       '
            s9 = '       '
            s10 = '       '
            readf, 1, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, format='(10A7)'  ;'(10D7.4)'
            if stregex(s1, '\*', /boolean) then s1 = 1.0
            if stregex(s2, '\*', /boolean) then s2 = 1.0
            if stregex(s3, '\*', /boolean) then s3 = 1.0
            if stregex(s4, '\*', /boolean) then s4 = 1.0
            if stregex(s5, '\*', /boolean) then s5 = 1.0
            if stregex(s6, '\*', /boolean) then s6 = 1.0
            if stregex(s7, '\*', /boolean) then s7 = 1.0
            if stregex(s8, '\*', /boolean) then s8 = 1.0
            if stregex(s9, '\*', /boolean) then s9 = 1.0
            if stregex(s10, '\*', /boolean) then s10 = 1.0
            spec[10*j:10*j+9,i] = float([s1, s2, s3, s4, s5, s6, s7, s8, s9, s10])
        endfor
        leftover = nlambda mod 10
        if leftover gt 0 then begin
            readf, 1, temp
            for k=0L,leftover-1 do begin
                s = '       '
                if k gt 0 then reads, temp, s, format='('+strcompress(7*k, /rem)+'X,D7.4)' else reads, temp, s, format='(A7)'
                
                if stregex(s, '\*', /boolean) then s = 1.0
                spec[10*j+k,i] = float(s)
            endfor
        endif
        i += 1
    endwhile
    abund = abund[0:i-1,*]
    spec = spec[*,0:i-1]
    w = where(spec gt 1.0, c)
    if c gt 0 then spec[w] = 1.0
    close, lun
    free_lun, lun
    spec = 1.0 - spec
    str = {elem:elem, abund:abund, lambda:lambda, spec:spec}
    if keyword_set(fits) then mwrfits, str, fitsfile, /create
    return, str
end
