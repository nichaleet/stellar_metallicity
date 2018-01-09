function make_spec_specify_template,grid_feh,grid_age,science
common sps_spec, sps, spsz, spsage
common get_sps, dlam, dataivar, datalam, wonfit, contmask, normalize,rest
	;setting
	redshift = science.zfit
	vdisp = science.vdisp
	reallambda = science.lambda 
	dlam = science.dlam
 
	znow = science.zfit
	won = where(science.fitmask eq 1 and finite(science.contdiv) and finite(science.contdivivar) and science.contdivivar gt 0 and reallambda/(1.+znow) gt 3500. and reallambda/(1.+znow) lt 7400., con)
	wonfit = won

	xmp = reallambda[won]
	ymp = science.contdiv[won]
	dymp = (science.contdivivar[won])^(-0.5)

	dataivar = science.telldivivar*(median(science.telldiv))^2
	datalam = reallambda
	contmask = science.contmask
	normalize = 0
	rest = 0
	
	arrsize= size(grid_feh)
	nfeh = arrsize[1]
	nage = arrsize[2]
	;loop over grid_feh and grid_age
	for iz=0,nfeh-1 do for ia=0,nage-1 do begin
		spsspec = get_sps_obs(xmp,[grid_feh[iz,ia],grid_age[iz,ia],vdisp,redshift])
		str = {spec:spsspec,feh:grid_feh[iz,ia],age:grid_age[iz,ia],vdisp:vdisp,redshift:redshift}		
		if iz eq 0 and ia eq 0 then strall=replicate(str,nfeh,nage)
		strall[iz,ia] = str
		if iz mod 50 eq 0 and ia mod 50 eq 0 then print, iz, ia
	endfor
	return,strall
end
