function transform_sdss_megacam, in, inerr, error=outerr, new=new, reverse=reverse
    a = 0.241
    b = 0.153
    c = 0.024
    d = keyword_set(new) ? 0.085 : 0.003
    e = 0.074

    matrix = [[1.0-a, a, 0.0, 0.0, 0.0], $
              [0.0, 1.0-b, b, 0.0, 0.0], $
              [0.0, -c, 1.0+c, 0.0, 0.0], $
              [0.0, 0.0, -d, 1.0+d, 0.0], $
              [0.0, 0.0, 0.0, -e, 1.0+e]]

    if keyword_set(reverse) then matrix = invert(matrix)
    out = matrix ## in

    outerr = dblarr(5)
    for i=0,4 do begin
        for j=0,4 do begin
            outerr[i] += (matrix[j,i]*inerr[j])^2.
        endfor
        outerr[i] = sqrt(outerr[i])
    endfor
    return, out
end
    
