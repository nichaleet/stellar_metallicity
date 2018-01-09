function make_spec_lowmass,grid_feh,grid_age
common sps_spec, sps, spsz, spsage
common get_sps, dlam, dataivar, datalam, wonfit, contmask, normalize,rest
	;setting
	redshift = 0.05
	vdisp = 150. ;FWHM in km/s 

	;get the obs setting
	;roughly for LRIS with R~1000
	;For 1" slit, the dispersion is 4 A. The detector is 0.135" per pix, so it's 0.54 A per pixel
	aperpix = 0.54
	resolution=1000.
	wlrange = [3700.,6200.]
	nlambda = fix((wlrange[1]-wlrange[0])/aperpix)
	datalam = findgen(nlambda)*aperpix+wlrange[0]
	dlam = datalam/resolution
		
	wonfit = where(datalam/(1.+redshift) gt 3500. and datalam/(1.+redshift) lt 7400., con)

	xmp = datalam[wonfit]

	normalize = 0
	rest = 0
	
	;loop over grid_feh and grid_age
	nfeh = n_elements(grid_feh)
	nage = n_elements(grid_age)
	for iz=0,nfeh-1 do for ia=0,nage-1 do begin
		spsspec = get_sps_obs(xmp,[grid_feh(iz),grid_age(ia),vdisp,redshift])
		str = {lambda:xmp,spec:spsspec,feh:grid_feh(iz),age:grid_age(ia),vdisp:vdisp,redshift:redshift}		
		if iz eq 0 and ia eq 0 then strall=replicate(str,nfeh,nage)
		strall[iz,ia] = str
		if iz mod 10 eq 0 and ia mod 10 eq 0 then print, iz, ia
	endfor
	return,strall
end
