pro combine_hires, imagenum, object, nowrite=nowrite
    n = n_elements(imagenum)

    for i=0,n-1 do begin
        numstr = string(imagenum[i], format=imagenum[i] lt 1000 ? '(I03)' : '(I04)')
        flux = mrdfits('Flux-'+numstr+'_1.fits', 0, hdr, /silent)
        err = mrdfits('Err-'+numstr+'_1.fits', 0, hdre, /silent)
        exptime = sxpar(hdr, 'EXPTIME')

        print, sxpar(hdr, 'OBJECT'), '   ', sxpar(hdr, 'DATE-OBS'), '    ', exptime

        fluxsum = i eq 0 ? flux : fluxsum + flux
        errsum = i eq 0 ? err^2d : errsum + err^2d
        exptimesum = i eq 0 ? exptime : exptimesum + exptime
    endfor
    errsum = sqrt(errsum)
    sxaddpar, hdr, 'EXPTIME', exptimesum
    sxaddpar, hdre, 'EXPTIME', exptimesum
    if ~keyword_set(nowrite) then mwrfits, fluxsum, object+'_1.fits', hdr, /create
    if ~keyword_set(nowrite) then mwrfits, errsum, object+'_1e.fits', hdre, /create

    for i=0,n-1 do begin
        numstr = string(imagenum[i], format=imagenum[i] lt 1000 ? '(I03)' : '(I04)')
        flux = mrdfits('Flux-'+numstr+'_2.fits', 0, hdr, /silent)
        err = mrdfits('Err-'+numstr+'_2.fits', 0, hdre, /silent)
        exptime = sxpar(hdr, 'EXPTIME')

        fluxsum = i eq 0 ? flux : fluxsum + flux
        errsum = i eq 0 ? err^2d : errsum + err^2d
        exptimesum = i eq 0 ? exptime : exptimesum + exptime
    endfor
    errsum = sqrt(errsum)
    sxaddpar, hdr, 'EXPTIME', exptimesum
    sxaddpar, hdre, 'EXPTIME', exptimesum
    if ~keyword_set(nowrite) then mwrfits, fluxsum, object+'_2.fits', hdr, /create
    if ~keyword_set(nowrite) then mwrfits, errsum, object+'_2e.fits', hdre, /create

    for i=0,n-1 do begin
        numstr = string(imagenum[i], format=imagenum[i] lt 1000 ? '(I03)' : '(I04)')
        flux = mrdfits('Flux-'+numstr+'_3.fits', 0, hdr, /silent)
        err = mrdfits('Err-'+numstr+'_3.fits', 0, hdre, /silent)
        exptime = sxpar(hdr, 'EXPTIME')

        fluxsum = i eq 0 ? flux : fluxsum + flux
        errsum = i eq 0 ? err^2d : errsum + err^2d
        exptimesum = i eq 0 ? exptime : exptimesum + exptime
    endfor
    errsum = sqrt(errsum)
    sxaddpar, hdr, 'EXPTIME', exptimesum
    sxaddpar, hdre, 'EXPTIME', exptimesum
    if ~keyword_set(nowrite) then mwrfits, fluxsum, object+'_3.fits', hdr, /create
    if ~keyword_set(nowrite) then mwrfits, errsum, object+'_3e.fits', hdre, /create
end
