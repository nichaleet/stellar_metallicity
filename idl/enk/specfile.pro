;+
; NAME:
;       specfile
;
; PURPOSE:
;       Return the spec1d filename(s) for a zresult structure
;       (array).
;
; CALLING SEQUENCE:
;       result = specfile(zcat)
;
; INPUT:
;       zcat      zresult structure, which can be an array
;
; RESULT:
;       string (array) of the spec1d filename
;
; OPTIONAL OUTPUTS:
;       found     indices for the zresult structure for which the
;                 spec1d files exist
;       notfound  indices for the zresult structure for which the
;                 spec1d files do not exist; note that these filenames
;                 still appear in the result
;-

function specfile, kz, found=found, notfound=notfound
    kztemp = kz
    obsdate = strcompress(kztemp.date)
    month_number = strmid(obsdate, 5, 2)
    month_word = strarr(n_elements(month_number))
    for i=0, n_elements(month_number)-1 do begin
        if month_number[i] eq '01' then month_word[i] = 'jan'
        if month_number[i] eq '02' then month_word[i] = 'feb'
        if month_number[i] eq '03' then month_word[i] = 'mar'
        if month_number[i] eq '04' then month_word[i] = 'apr'
        if month_number[i] eq '05' then month_word[i] = 'may'
        if month_number[i] eq '06' then month_word[i] = 'jun'
        if month_number[i] eq '07' then month_word[i] = 'jul'
        if month_number[i] eq '08' then month_word[i] = 'aug'
        if month_number[i] eq '09' then month_word[i] = 'sep'
        if month_number[i] eq '10' then month_word[i] = 'oct'
        if month_number[i] eq '11' then month_word[i] = 'nov'
        if month_number[i] eq '12' then month_word[i] = 'dec'        
    endfor
    date = strmid(obsdate, 0, 4) + month_word + strmid(obsdate, 8, 2)
    gooddatew = where(strlen(date) eq 9, cw)
    if cw ne n_elements(kztemp) then message, 'One or more dates are invalid.'
    date = date[gooddatew]
    kztemp = kztemp[gooddatew]

    specfile = getenv('D2_RESULTS') + '/' + kztemp.maskname + '/' + date $
      + '/spec1d.' + kztemp.maskname + '.' + $
      string(long(kztemp.slitname), format='(I03)') + $
      '.' + string(long(kztemp.objname), format='(I8)') + '.fits'
    notfound = [-1]
    for i=0,n_elements(specfile)-1 do begin
        if not(file_test(specfile[i])) then begin
            print, specfile[i] + ' not found'
            notfound = [notfound, i]
        endif
    endfor
    n = n_elements(notfound)
    if n gt 1 then notfound = notfound[1:n-1]
    found = complement(notfound, n_elements(specfile))
    return, specfile
end
