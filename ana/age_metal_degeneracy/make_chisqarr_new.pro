function make_chisqarr_new,fehref,ageref,sn
common degen, spec_Arr,grid_feh,grid_age

	;find where is the spec closest to the input age and feh
        min_dist = min(sqrt((spec_arr.feh-fehref)^2+(spec_arr.age-ageref)^2),iref)
	if min_dist gt 0.05 then stop,'off grid'

	;make noisy ref spectra
	npix = n_elements(spec_arr(iref).lambda)
	refspec_err = abs(spec_arr(iref).spec/sn)
	refspec = spec_arr(iref).spec+randomn(seed,npix)*refspec_err

	;make chisqarr
	chisqarr = fltarr(size(grid_feh,/dimensions))
	for i=0,n_elements(spec_arr)-1 do begin
		loc = where(grid_feh eq spec_arr[i].feh and grid_age eq spec_arr[i].age,cloc)
		if cloc ne 1 then stop,'oops'
                ind = array_indices(grid_feh,loc)
		chisqarr[ind[0],ind[1]] = total((refspec-spec_arr[i].spec)^2/refspec_err)	
	endfor
	
	return,chisqarr
end
