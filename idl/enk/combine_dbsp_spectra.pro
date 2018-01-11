pro combine_dbsp_spectra, bluenums, rednums, filename, overwrite=overwrite, aper=aper
    if ~keyword_set(aper) then aper = 1
    filename += '.fits'
    nblue = n_elements(bluenums)
    nred = n_elements(rednums)
    if file_test(filename) and ~keyword_set(overwrite) then message, filename+' exists.'
    for i=0,nblue-1 do begin
        blue = mrdfits('blue'+string(bluenums[i], format='(I04)')+'b.ms.fits', 0, hdrb, /silent)
        lambdablue = (dindgen(n_elements(blue[*,0,0]))+1.0-double(sxpar(hdrb, 'CRPIX1')))*double(sxpar(hdrb, 'CD1_1')) + double(sxpar(hdrb, 'CRVAL1'))
        fluxblue = blue[*,aper-1,0]
        sigmablue = blue[*,aper-1,3]
        
        if i eq 0 then begin
            lambdab = lambdablue
            fluxb = fluxblue
            sigmab = sigmablue^2.
        endif else begin
            x_specrebin, lambdablue, fluxblue, lambdab, fluxbtemp, /silent
            fluxb += fluxbtemp
            x_specrebin, lambdablue, sigmablue, lambdab, sigmabtemp, /silent
            sigmab += sigmabtemp^2.
        endelse
    endfor
    for i=0,nred-1 do begin
        red = mrdfits('red'+string(rednums[i], format='(I04)')+'b.ms.fits', 0, hdrr, /silent)
        lambdared = (dindgen(n_elements(red[*,0,0]))+1.0-double(sxpar(hdrr, 'CRPIX1')))*double(sxpar(hdrr, 'CD1_1')) + double(sxpar(hdrr, 'CRVAL1'))
        fluxred = red[*,aper-1,0]
        sigmared = red[*,aper-1,3]
        
        if i eq 0 then begin
            lambdar = lambdared
            fluxr = fluxred
            sigmar = sigmared^2.
        endif else begin
            x_specrebin, lambdared, fluxred, lambdar, fluxrtemp, /silent
            fluxr += fluxrtemp
            x_specrebin, lambdared, sigmared, lambdar, sigmartemp, /silent
            sigmar += sigmartemp^2.
        endelse
    endfor
    sigmab = sqrt(sigmab)
    sigmar = sqrt(sigmar)

    spec = {lambda:[lambdab, lambdar], spec:[fluxb, fluxr], err:[sigmab, sigmar]}
    mwrfits, spec, filename, /create
end
