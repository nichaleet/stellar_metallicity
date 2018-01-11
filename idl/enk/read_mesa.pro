function read_mesa, filename
    varnames = ' '
    openr, lun, filename, /get_lun
    skip_lun, lun, 5, /lines
    readf, lun, varnames, format='(A3000)'
    varnames = strsplit(varnames, /extract)
    nvars = n_elements(varnames)
    
    str0s = replicate(', 0d', nvars)
    str0s = strjoin(str0s)
    exresult = execute('data = create_struct(varnames'+str0s+')')
    data = replicate(data, 100000)

    strin = ' '
    i = 0L
    while ~eof(lun) do begin
        readf, lun, strin, format='(A3000)'
        vars = strsplit(strin, /extract)
        for j=0,nvars-1 do begin
            data[i].(j) = double(vars[j])
        endfor
        i++
    endwhile

    close, lun
    free_lun, lun

    data = data[0:i-1]
    return, data
end
