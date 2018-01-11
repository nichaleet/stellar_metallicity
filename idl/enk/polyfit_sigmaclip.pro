function polyfit_sigmaclip, x, y, yerr, degree, sigclip=sigclip, _extra=extra
    good = 0
    if ~keyword_set(sigclip) then sigclip = 3.0
    w = indgen(n_elements(y))
    cprev = -1
    while good eq 0 do begin
        p = poly_fit(x[w], y[w], degree, measure_errors=yerr[w], yfit=yfit, extra)
        w = where(abs(y[w]-yfit) lt sigclip*yerr[w], c)
        if c lt degree then begin
            message, 'Not enough points left for a quadratic fit.  No sigma clipping performed.', /info
            p = poly_fit(x, y, degree, measure_errors=yerr, yfit=yfit, extra)
            good = 1
        endif
        if c eq cprev then good = 1
        cprev = c
    endwhile
    return, p
end
