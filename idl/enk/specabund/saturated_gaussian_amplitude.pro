function saturated_gaussian_amplitude, ew_in, fwhm
    common satgauss, a, ew
    ;sigma = fwhm/(2d * sqrt(2d * alog(2d + a)))
    sigma = abs(fwhm) / 2.354820045
    ewsigma = ew*sigma
    w = where(ewsigma ge ew_in-0.01 and ewsigma le ew_in+0.01, c)
    ;print, ew_in, c
    amplitude = interpolate(a[w], findex(ewsigma[w], ew_in))
    return, amplitude[0]
end
