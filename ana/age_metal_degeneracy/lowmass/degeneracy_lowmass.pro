pro degeneracy_lowmass,sn=sn,redo=redo
common degen, spec_Arr, grid_feh, grid_age
	if n_Elements(redo) eq 1 then begin
		;make spectra with these grids
		grid_feh = findgen(126)/75.-1.5 ;range from [-1.5,0.17] dex 
		grid_age = findgen(121)/15.+0.1 ;range from [0.1,8.1] Gyr
		spec_arr = make_spec_lowmass(grid_feh,grid_age)
		mwrfits,spec_arr,'grid_spec_lris.fits',/silent,/create
	endif else spec_arr = mrdfits('grid_spec_lris.fits',1)
	fehpos = uniq(spec_arr.feh,sort(Spec_arr.feh))
	grid_feh = spec_Arr(fehpos).feh
	agepos = uniq(spec_arr.age,sort(Spec_arr.age))
	grid_age = spec_Arr(agepos).age

	;find chisquare at each reference point
	ref_ia = [20,30,40,50,60,70,80,90,100]+1
	ref_iz = [20,30,40,50,60,70,80,90,100,110]+1 
	for i=0,n_Elements(ref_ia)-1 do for j=0,n_elements(ref_iz)-1 do begin
		ageref = grid_age(ref_ia(i))
		fehref = grid_feh(ref_iz(j))
		uncert = make_chisqarr(ref_ia(i),ref_iz(j),sn)
	endfor
	;stop
end
