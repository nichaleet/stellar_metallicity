function saturated_gaussian, x, a
    fwhm = a[1]
    dlambda = a[2]
    amplitude = saturated_gaussian_amplitude(1d-3*a[0], fwhm)
    sigma = fwhm / (2d * sqrt(2d * alog(2d + 0*amplitude)))
    ;g = amplitude*exp(-(x^2d / (2d * sigma^2d)))
    g = amplitude*exp(-((x-dlambda)^2d / (2d * sigma^2d)))
    func = g/(1d + g/0.95)
    func += a[3]
    return, func
end
