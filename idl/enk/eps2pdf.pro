pro eps2pdf
    files = file_search('*.eps', count=c)
    if c gt 0 then begin
        for i=0,c-1 do begin
            print, strcompress(files[i], /rem)
            spawn, 'ps2pdf '+strcompress(files[i], /rem)
        endfor
    endif
end
