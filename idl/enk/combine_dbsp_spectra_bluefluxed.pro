pro combine_dbsp_spectra_bluefluxed, bluenums, basename, overwrite=overwrite, aper=aper
    if ~keyword_set(aper) then aper = 1
    filename = basename+'_fluxed.fits'
    nblue = n_elements(bluenums)
    nred = n_elements(rednums)
    if file_test(filename) and ~keyword_set(overwrite) then message, filename+' exists.'

    counts = mrdfits(basename+'_blue.fits', 0, hdrb, /silent)
    lambda1 = (dindgen(2834)+1.0-double(sxpar(hdrb, 'CRPIX1')))*double(sxpar(hdrb, 'CD1_1')) + double(sxpar(hdrb, 'CRVAL1'))
    fluxed = mrdfits(basename+'_blue_fluxed.fits', 0, hdrb, /silent)
    lambda2 = (dindgen(2834)+1.0-double(sxpar(hdrb, 'CRPIX1')))*double(sxpar(hdrb, 'CD1_1')) + double(sxpar(hdrb, 'CRVAL1'))
    if ~array_equal(lambda1, lambda2) then message, "There's a problem here."

    for i=0,nblue-1 do begin
        blue = mrdfits('blue'+string(bluenums[i], format='(I04)')+'b.ms.fits', 0, hdrb, /silent)
        lambdablue = (dindgen(2834)+1.0-double(sxpar(hdrb, 'CRPIX1')))*double(sxpar(hdrb, 'CD1_1')) + double(sxpar(hdrb, 'CRVAL1'))
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
    sigmab = sqrt(sigmab)

    x_specrebin, lambda1, fluxed/counts, lambdab, fluxfactor, /silent
    fluxb *= fluxfactor
    sigmab *= fluxfactor

    spec = {lambda:lambdab, spec:fluxb, err:sigmab}
    mwrfits, spec, filename, /create
end
