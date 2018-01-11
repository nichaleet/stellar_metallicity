function sigma_clip, array, sigma, iter
    ;if ~arg_present(sigma) then sigma = 2
    ;if ~arg_present(iter) then iter = 5
    n = n_elements(array)
    w = lindgen(n)
    for i=0,iter-1 do begin
        if n lt 2 then begin
            message, 'The array contains fewer than 2 elements.', /info
            break
        endif
        mean = mean(array[w])
        sd = stddev(array[w])
        wnew = w[where(abs(array[w]-mean) le sd, n)]
        if array_equal(wnew, w) then break else w = wnew
    endfor
    return, w
end
