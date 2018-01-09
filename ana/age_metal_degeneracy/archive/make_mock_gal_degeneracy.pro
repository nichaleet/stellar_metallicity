pro make_mock_gal_degeneracy
;This program is a nest inside test_degeneracy_formation_time.pro
;This is the part where it makes a fits files of mock galaxies whose the real MZR follows
;a tight MZR but scattered by age-metal degeneracy

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;setting
	massrange = [9.0,11.5]
	ngals = 200.
	SNrange = [15,25] ;SDSS samples
	z=0.05
	mzrparam=[9.01704,9.91704,0.817038]
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	indi = {sn:0.,mass:0.,fehideal:0.,ageideal:0.,feh:0.,feherr:0.,age:0.,ageerr:0.}
	str  = replicate(indi,ngals)
	
	;assign mass randomly within the range
	str.mass=randomu(seed,ngals)*(massrange[1]-massrange[0])+massrange[0]
	;assign ideal metallicities
	str.fehideal = mzrparam[0]-alog10(1.+10.^((str.mass-mzrparam[1])*(-1.*mzrparam[2])))-8.9
	;get ideal age randomly within a range
	for i=0,ngals-1 do begin
		massnow = str(i).mass
		minage = 1. ;1 Gyr
		maxage = 5.*massnow-41.< 10 ;function of mass but less than 10 Gyr
		;this is from fitting linear fit by eye the maximum range of age as a function of mass from age-mass plot
		str(i).ageideal = randomu(seed)*(maxage-minage)+minage
	endfor
	;assign random sn
	str.sn=randomu(seed,ngals)*(snrange[1]-snrange[0])+snrange[0]
		
	;read the synthesis spec of sdss
	spec_arr = mrdfits('new_grid_spec_sdss.fits',1)
	;get values of age metallciities
	fehpos = uniq(spec_arr.feh,sort(Spec_arr.feh))
	grid_feh = spec_Arr(fehpos).feh
	agepos = uniq(spec_arr.age,sort(Spec_arr.age))
	grid_age = spec_Arr(agepos).age
	nfeh = n_Elements(grid_feh)
	nage = n_Elements(grid_age)
	;make 2d arrays from 1d arrays
	grid_Feh = rebin(grid_feh,nfeh,nage)
	grid_age = transpose(rebin(grid_age,nage,nfeh)) ;range from [0.5,13] Gyr

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;run a loop for each 'galaxy' to get the 'observed feh and observed age'
	for i=0,ngals-1 do begin
		;find the matched spectrum
		dist = abs(spec_arr.feh-str(i).fehideal)+abs(spec_arr.age-str(i).ageideal)
		loc = where(dist eq min(dist))
		print,strtrim(string(fix(i)),2)+' FeH=',str(i).fehideal,' use ',spec_arr(loc).feh
		print,strtrim(string(fix(i)),2)+' AGE=',str(i).ageideal,' use ',spec_arr(loc).age
		;create noised spectrum
		refspecstr = spec_arr(loc)
		npix = n_elements(refspecstr.lambda)
		refspec_err = abs(refspecstr.spec/str(i).sn)
		refspec = refspecstr.spec+randomn(seed,npix)*refspec_err
		;make chisq arr
		chisqarr = fltarr(100,126)           
		for j=0,n_elements(spec_arr)-1 do begin
			loc = where(grid_feh eq spec_arr[j].feh and grid_age eq spec_arr[j].age,cloc)
			if cloc ne 1 then stop,'you have a problem'
			chisqarr(loc) = total((refspec-spec_arr[j].spec)^2/refspec_err)
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
	mwrfits,str,'mock_gal_degeneracy.fits',/silent,/create
end


end
