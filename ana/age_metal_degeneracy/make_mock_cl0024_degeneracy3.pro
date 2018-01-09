pro make_mock_cl0024_degeneracy3
;This program is a nest inside test_degeneracy_formation_time.pro
;This is the part where it makes a fits files of mock galaxies whose the real MZR follows
;a tight MZR but scattered by age-metal degeneracy
common mzr_lin_param, mzrparam
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;setting
	sci = mrdfits('/scr2/nichal/workspace2/sps_fit/data/all_cl0024/sps_fit.fits.gz',1,/silent)	
	cat = mrdfits('/scr2/nichal/keck/deimos/Cl0024MS0451/Cl0024master.v7.fits.gz',1,/silent)

	file_matchedcat = '/scr2/nichal/workspace2/ana/cl0024/cl0024_matchedcat.sav'
	restore,file_matchedcat
	cat = cat(matched_cat)
	nsci = n_elements(sci)
	;check the match
	for i=0,nsci-1 do begin
	        if randomu(seed) lt 0.2 then begin
	                gcirc,2,cat(i).ra,cat(i).dec,sci(i).ra,sci(i).dec,dis
	                if dis gt 1 then stop,'check the match of ra/dec'  ;if the distance is greater than 1 arcsecond
	        endif
	endfor
	print,'Roughly checked the match between the catalog and science catalog.'

	wgoodmem = where(sci.goodfit+sci.good eq 2 and sci.zfit ge 0.38 and sci.zfit le 4.0 and sci.logmstar gt 7., cgals)
	wbadmem  = where(sci.goodfit+sci.good eq 1 and sci.zfit ge 0.38 and sci.zfit le 4.0 and sci.feh ne -999. and sci.logmstar gt 7., cbadgals)
	wallmem = [wgoodmem,wbadmem]
	wearlymem = where(cat(wallmem).morph le 2 and cat(wallmem).morph ge 0,cearlymem)
	wearlymem = wallmem(wearlymem)	
	sci = sci(wearlymem)

	ngals = n_elements(sci)
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	indi = {sn:0,mass:0.,fehideal:0.,ageideal:0.,feh:0.,feherr:0.,age:0.,ageerr:0.,z:0.,vdisp:0.}
	str  = replicate(indi,ngals)

	;assign mass as the observed mass
	str.mass = sci.logmstar
	;assign ideal metallicities
	str.fehideal = mzrparam[0]+str.mass*mzrparam[1]
	;get ideal age as the observed age
	str.ageideal = sci.age
	str.z = sci.zfit
	str.sn = sci.sn
	str.vdisp = sci.vdisp

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;setting
	nfeh = 120
	nage = 126
	grid_Feh = rebin(findgen(nfeh)/120.-0.8,nfeh,nage) ;range from [-0.8,0.192] dex ;
	grid_Age = transpose(rebin(findgen(126)/10.+0.5 ,nage,nfeh)) ;range from [0.5,13] Gyr
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;run a loop for each 'galaxy' to get the 'observed feh and observed age'
	for i=0,ngals-1 do begin
		print, 'now doing cl0024 galaxy number', i
		;create synthesis spec with template from sci
		spec_arr = make_spec_specify_template(grid_feh,grid_age,sci(i))          

		;find the matched spectrum
		dist = abs(spec_arr.feh-str(i).fehideal)+abs(spec_arr.age-str(i).ageideal)
		loc = where(dist eq min(dist))
		print,strtrim(string(fix(i)),2)+' FeH=',str(i).fehideal,' use ',spec_arr(loc).feh
		print,strtrim(string(fix(i)),2)+' AGE=',str(i).ageideal,' use ',spec_arr(loc).age

		;create noised spectrum
		;use the noise as a fn of wl from observed spectra
		refspecstr = spec_arr(loc)
		won = where(sci(i).fitmask eq 1 and finite(sci(i).contdiv) and finite(sci(i).contdivivar) and sci(i).contdivivar gt 0 and sci(i).lambda/(1.+sci(i).zfit) gt 3500. and sci(i).lambda/(1.+sci(i).zfit) lt 7400., con)
		if n_Elements(refspecstr.spec) ne con then stop,'oh oh'
		refspec_err = abs(refspecstr.spec/(sci(i).contdiv(won)*Sqrt(sci(i).contdivivar(won))))
		npix=con
		refspec = refspecstr.spec+randomn(seed,npix)*refspec_err

		;make chisq arr
		chisqarr = fltarr(nfeh,nage)           
		for iz=0,nfeh-1 do for ia=0,nage-1 do begin
			if grid_feh[iz,ia] ne spec_arr[iz,ia].feh or grid_age[iz,ia] ne spec_arr[iz,ia].age then stop,'emm sth is wrong'
			chisqarr[iz,ia] = total((refspec-spec_arr[iz,ia].spec)^2/refspec_err)
		endfor	
		;select the points within 2 sigmas (chisq<=minchisq+4)
		wgoodchi = where(chisqarr lt min(chisqarr)+4.,cgoodchi)
		goodchi = chisqarr(wgoodchi)
		goodchi = goodchi-max(goodchi)
		if cgoodchi lt 3 then stop,'eh, something is wrong'
		;get the probabilities of these guys
		prob = exp(-0.5*goodchi)		
		cumuprob =total(prob,/cumulative)
		randomnum = randomu(seed)*cumuprob(cgoodchi-1)
		chosen = value_locate(cumuprob,randomnum)+1
		str(i).feh = grid_feh(wgoodchi(chosen))
		str(i).age = grid_age(wgoodchi(chosen))
		;now do the sigmas
		wgoodchi = where(chisqarr lt min(chisqarr)+1.,cgoodchi)
		str(i).feherr = 0.5*(max(grid_Feh(wgoodchi))-min(grid_Feh(wgoodchi)))	
		str(i).ageerr = 0.5*(max(grid_age(wgoodchi))-min(grid_age(wgoodchi)))
	endfor
	mwrfits,str,'mock_cl0024_degeneracy3.fits',/silent,/create
end


