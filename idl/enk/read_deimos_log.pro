pro deimoslog__define
    deimoslog = {deimoslog, date:' ', filename:' ', objname:' ', pa:0d, lamps:' ', mask:' ', grating:' ', cwave:0d, filter:' ', focus:0d, exptime:0d, airmass:0d, comment:' ', comment2:' ', seeing:0d}
end


function read_deimos_log, date
    logfile = strmid(date, 3, 4, /reverse_offset) ne '.log' ? getenv('DEIMOS_DATA')+'/logfiles/'+date+'.log' : date
    date = file_basename(logfile, '.log')
    deimoslog = replicate({deimoslog}, 1000)
    n = 0L
    openr, lun, logfile, /get_lun
    skip_lun, lun, 7, /lines
    instring = ' '
    comment2 = ''
    seeing = 0.0

    while ~eof(lun) do begin
        readf, lun, instring, format='(A)'

        if strmid(instring, 0, 1) eq '#' then begin
            comment2 = instring
            if stregex(instring, 'seeing', /boolean) then begin
                s1 = strpos(instring, '=')+1
                s2 = strpos(instring, '"')
                if s2 ge 0 then begin
                    if s1 le 1 then s1 = strpos(instring, '.', s2+1, /reverse_search)-1
                    if s1 le 1 then s1 = strpos(instring, '~', s2+1, /reverse_search)+1
                    if s1 le 1 then s1 = strpos(instring, '>', s2+1, /reverse_search)+1
                    if s1 le 1 then s1 = strpos(instring, '<', s2+1, /reverse_search)+1
                    if s1 le 1 then s1 = strpos(instring, ' ', s2+1, /reverse_search)+1
                    s1 >= 0
                    s1 >= s2-6
                    seeing = double(strmid(instring, s1, s2-s1))
                endif
            endif
        endif
        if strmid(instring, 0, 1) ne 'd' then continue
        if strmid(instring, 10, 1) ne ' ' then continue

        ;print, instring
        deimoslog[n].date = date
        deimoslog[n].filename = strtrim(strmid(instring, 0, 10), 2)
        deimoslog[n].objname = strtrim(strmid(instring, 11, 15), 2)
        pastring = strmid(instring, 26, 7)
        deimoslog[n].pa = strtrim(pastring, 2) eq '?' ? -999.0 : double(pastring)
        stringrempos = strsplit(strmid(instring, 34))
        stringrem = strsplit(strmid(instring, 34), /extract)
        if strmid(instring, 34, 1) eq ' ' then begin
            stringrempos = [0, stringrempos]
            stringrem = [' ', stringrem]
        endif
        ns = n_elements(stringrem)
        deimoslog[n].lamps = strtrim(stringrem[0], 2)
        deimoslog[n].mask = strtrim(stringrem[1], 2)
        deimoslog[n].grating = strtrim(stringrem[2], 2)
        if stringrempos[3]-stringrempos[2] gt 14 then begin
            deimoslog[n].cwave = 0.0
            e = 0
        endif else begin
            deimoslog[n].cwave = double(strtrim(stringrem[3], 2))
            e = 1
        endelse
        deimoslog[n].filter = strtrim(stringrem[3+e], 2)
        focusstring = strtrim(stringrem[4+e], 2)
        deimoslog[n].focus = focusstring eq '?' ? -999.0 : double(focusstring)
        exptimestring = strtrim(stringrem[5+e], 2)
        deimoslog[n].exptime = exptimestring eq '?' ? -999.0 : double(exptimestring)
        airmassstring = strtrim(stringrem[6+e], 2)
        deimoslog[n].airmass = airmassstring eq '?' ? -999.0 : double(airmassstring)
        if ns gt 7+e then deimoslog[n].comment = strtrim(strjoin(stringrem[7+e:ns-1], ' '), 2)

        deimoslog[n].comment2 = comment2
        comment2 = ''
        deimoslog[n].seeing = seeing

        n++
    endwhile
    deimoslog = deimoslog[0:n-1]

    free_lun, lun
    close, lun

    return, deimoslog
end
