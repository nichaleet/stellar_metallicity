function gauss_fixed_area, x, a, ew=ew
    amplitude = 1d-3*ew / (a[0]*sqrt(!DPI))
    f = amplitude*exp(-1d*(x/a[0])^2d)
    return, f
end


pro plot_daospec, star, directory=directory, overplot=overplot, ew=ew, vr=vr
    if ~keyword_set(directory) then directory = '.'
    contfile = file_search(directory+'/'+star+'*C.fits', count=c)
    for i=0,c-1 do begin
        specfile = strmid(contfile[i], 0, strlen(contfile[i])-6)+'.fits'
        residfile = strmid(contfile[i], 0, strlen(contfile[i])-6)+'R.fits'
        speci = mrdfits(specfile, 0, hdrs, /silent)
        conti = mrdfits(contfile[i], 0, hdrc, /silent)
        residi = mrdfits(residfile, 0, hdrr, /silent)
        npix = sxpar(hdrs, 'NAXIS1')
        dlambda = sxpar(hdrs, 'CDELT1')
        lambda0 = sxpar(hdrs, 'CRVAL1')
        lambdai = dindgen(npix)*dlambda+lambda0
        
        if i eq 0 then begin
            lambda = lambdai
            spec = speci / conti
            resid = residi
        endif else begin
            lambda = [lambda, lambdai]
            spec = [spec, speci / conti]
            resid = [resid, residi]
        endelse
    endfor

    if ~keyword_set(vr) then begin
        restore, getenv('CALTECH')+'hires/n5024/'+star+'_vr.sav'
    endif
    c = 2.99792458d5
    lambda /= 1d + vr/c
    
    vsym, 24
    if keyword_set(overplot) then soplot, lambda, spec, yrange=[0, 2], symsize=0.4 else splot, lambda, spec, psym=-8, yrange=[0, 2], symsize=0.4
    ;if keyword_set(overplot) then soplot, lambda, resid+1d, color=fsc_color('red'), symsize=0.4 else soplot, lambda, resid+1d, psym=8, color=fsc_color('red'), symsize=0.4

    if keyword_set(ew) then begin
        readcol, getenv('CALTECH')+'hires/n5024/'+star+'/'+star+'.ew', lambdaw, speciesw, epw, loggfw, dampingw, eww, format='D,D,D,D,D,D'
        ews = {lambda:0d, species:0d, ep:0d, loggf:0d, ew:0d}
        ews = replicate(ews, n_elements(lambdaw))
        ews.lambda = lambdaw
        ews.species = speciesw
        ews.ep = epw
        ews.loggf = loggfw
        ews.ew = eww

        halfwindow = 0.25
        for i=0,n_elements(lambdaw)-1 do begin
            w = where(lambda gt ews[i].lambda-halfwindow and lambda lt ews[i].lambda+halfwindow, c)
            if c lt 10 then continue

            a = mpfitfun('gauss_fixed_area', lambda[w]-ews[i].lambda, 1.0-spec[w], replicate(1d, c), [6.0*0.03], /quiet, functargs={ew:ews[i].ew})
            sigma = a[0]
            ;sigma = 7.0*0.03/2.235*ews[i].lambda/5170d
            amplitude = 1d-3*ews[i].ew / (sigma*sqrt(!DPI))
            g = amplitude*exp(-1d*((lambda[w]-ews[i].lambda)/sigma)^2d)
            ;print, sigma

            soplot, lambda[w], 1.0-g, color=fsc_color('red'), thick=3
        endfor
    endif
end
