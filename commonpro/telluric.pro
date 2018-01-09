pro telluric,str
    tell = mrdfits('/scr2/nichal/workspace2/telluric/deimos_telluric_1.0.fits', 1, /silent)
    wtell = n_elements(tell)-1 ;the last file
    tell = tell[wtell]
    aratio = str.airmass/((tell).airmass)
    telllambda = (tell)[0].lambda
    tellspec = ((tell)[0].spec)^aratio
    tellivar = ((tell)[0].ivar)*(((tell)[0].spec)/(aratio*tellspec))^2.
    ivarmissing = 10d10
    w = where(tellivar ge 1d8, c) ;these are not in the telluric band
    if c gt 0 then tellivar[w] = ivarmissing
    f = where(finite(tellspec) and tellspec gt 0 and finite(tellivar) and tellivar gt 0 and tellivar lt 1d8) ;where the band is
    telllambda = telllambda[f]
    tellspec = tellspec[f]
    tellivarfull = tellivar
    tellivar = tellivar[f]

    tellspecnew = interpolate(tellspec, findex(telllambda, str.lambda), missing=1.)
    tellivarnew = interpolate(tellivarfull, findex(telllambda, str.lambda), missing=ivarmissing)
    
    npix = n_elements(str.lambda)
    str=create_struct(str,'telldiv',dblarr(npix))
    str=create_struct(str,'telldivivar',dblarr(npix))
    str.telldiv = str.spec / tellspecnew
    str.telldivivar = (str.ivar^(-1.) + (str.telldiv)^2.*tellivarnew^(-1.))^(-1.) * tellspecnew^2.
end
