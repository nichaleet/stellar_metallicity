pro abund__define
    abund = {abund, element:' ', ion:0, lambda:0., ep:0., loggf:0., ew:0., rw:0., abund:0., upperlimit:0}
end


function read_moog, outfile, epslope=epslope, epint=epint, epcorr=epcorr, rwslope=rwslope, rwint=rwint, rwcorr=rwcorr, lambdaslope=lambdaslope, lambdaint=lambdaint, lambdacorr=lambdacorr, upperlimit=upperlimit
    openr, lun, outfile, /get_lun
    instring = ' '
    readabund = 0
    abund = replicate({abund}, 1000)
    if keyword_set(upperlimit) then abund.upperlimit = 1
    i = 0L
    while ~eof(lun) do begin
        readf, lun, instring, format='(A80)'
        case strmid(instring, 0, 4) of
            'Abun': begin
                species = strmid(instring, 30, 5)
                readabund = 0
            end
            'wave': readabund = 1
            'E.P.': begin
                if species eq 'Fe I ' then begin
                    epslope = double(strmid(instring, 27, 7))
                    epint = double(strmid(instring, 48, 7))
                    epcorr = double(strmid(instring, 72, 7))
                endif
                readabund = 0
            end
            'R.W.': begin
                if species eq 'Fe I ' then begin
                    rwslope = double(strmid(instring, 27, 7))
                    rwint = double(strmid(instring, 48, 7))
                    rwcorr = double(strmid(instring, 72, 7))
                endif
                readabund = 0
            end
            'wav.': begin
                if species eq 'Fe I ' then begin
                    lambdaslope = double(strmid(instring, 27, 7))
                    lambdaint = double(strmid(instring, 48, 7))
                    lambdacorr = double(strmid(instring, 72, 7))
                endif
                readabund = 0
            end
            'aver': readabund = 0
            ' No ': readabund = 0
            else: begin
                if readabund then begin
                    element = strtrim(strmid(species, 0, 2), 2)
                    ionstr = strtrim(strmid(species, 3, 2), 2)
                    case ionstr of
                        'I': ion = 1
                        'II': ion = 2
                        'III': ion = 3
                    endcase
                    abund[i].element = element
                    abund[i].ion = ion
                    abund[i].lambda = double(strmid(instring, 3, 7))
                    abund[i].ep = double(strmid(instring, 13, 7))
                    abund[i].loggf = double(strmid(instring, 23, 7))
                    abund[i].ew = double(strmid(instring, 33, 7))
                    abund[i].rw = double(strmid(instring, 43, 7))
                    abund[i].abund = double(strmid(instring, 53, 7))
                    i++
                endif else readabund = 0
            end
        endcase
    endwhile
    close, lun
    free_lun, lun
    abund = abund[0:i-1]
    ;abund = abund[where(abund.abund lt 100, c)]
    c = n_elements(abund)

    for j=c-1,0,-1 do begin
        w = where(abund.element eq abund[j].element and abund.lambda eq abund[j].lambda and abund.ep eq abund[j].ep and abund.loggf eq abund[j].loggf, cw)
        if cw gt 1 then begin
            w = w[where(w ne j)]
            abund = abund[complement(w, n_elements(abund))]
            j -= 1
        endif
    endfor
    return, abund
end
