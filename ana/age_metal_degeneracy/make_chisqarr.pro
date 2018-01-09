function make_chisqarr,refa,refz,sn
common degen, spec_Arr, grid_feh, grid_age
    ageref = grid_age(refa)
    fehref = grid_feh(refz)

	;find where is the ref spec
	iref = where(spec_arr.feh eq fehref and spec_arr.age eq ageref,ciref)
	if ciref ne 1 then stop
	refspecstr = spec_arr(iref)
	;make noisy ref spectra
	npix = n_elements(refspecstr.lambda)
	refspec_err = abs(refspecstr.spec/sn)
	refspec = refspecstr.spec+randomn(seed,npix)*refspec_err

	;make chisqarr
	chisqarr = findgen(n_elements(grid_feh),n_elements(grid_age))
	feharr = chisqarr
	agearr = chisqarr

	for i=0,n_elements(spec_arr)-1 do begin
		iz = where(grid_feh eq spec_arr[i].feh,ciz)
		ia = where(grid_age eq spec_arr[i].age,cia)
		if ciz ne 1 or cia ne 1 then stop
		chisqarr(iz,ia) = total((refspec-spec_arr[i].spec)^2/refspec_err)				
		feharr(iz,ia) = grid_feh(iz)
		agearr(iz,ia) = grid_age(ia)
	endfor
	;find min max of uncertainties
	minchi = min(chisqarr)
	inrange = where(chisqarr lt minchi+1.,cinrange)
	if cinrange lt 2 then stop
	ageuncertain = minmax(agearr(inrange))
	fehuncertain = minmax(feharr(inrange))	

	;write fits file
	mkhdr,header,chisqarr
	paraname = ['CTYPE1','CTYPE2']
	paraval =['[Fe/H]','Age']
	for ii=0,1 do sxaddpar,header,paraname(ii),paraval(ii)
	paraname = ['CRVAL1','CDELT1','CRPIX1','CROTA1','CRVAL2','CDELT2','CRPIX2','CROTA2']
	paraval = [grid_feh(0),grid_feh(1)-grid_feh(0),0.,0.,grid_age(0),grid_age(1)-grid_age(0),0.,0.]
	for ii =0, 7 do sxaddpar,header,paraname(ii),paraval(ii)
	dir='sdssoutputsn'+strtrim(string(fix(sn)),2)
	if file_test(dir,/directory) eq 0 then file_mkdir,dir
	writefits,dir+'/chisq_age'+strtrim(string(fix(ageref*10)),2)+'_feh'+strtrim(string(fix(fehref*100)),2)+'.fits',chisqarr,header
	
	return,{agerange:ageuncertain,fehrange:fehuncertain}
end
