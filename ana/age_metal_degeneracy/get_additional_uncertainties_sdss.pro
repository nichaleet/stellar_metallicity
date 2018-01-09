pro get_additional_uncertainties_sdss,feh,age,sn,objname,choi=choi
common degen, spec_arr,grid_feh,grid_age
   spec_arr = mrdfits('/scr2/nichal/workspace2/ana/age_metal_degeneracy/new_grid_spec_sdss.fits',1)
   feharr = spec_arr(uniq(spec_arr.feh,sort(spec_arr.feh))).feh 
   agearr = spec_arr(uniq(spec_arr.age,sort(spec_arr.age))).age
   nfeh = n_elements(feharr)
   nage = n_elements(agearr)
   grid_feh = rebin(feharr,nfeh,nage)
   grid_age = transpose(rebin(agearr,nage,nfeh))

   ngals = n_elements(feh)
   chisq_cube = fltarr(nfeh,nage,ngals)
   ageupper = fltarr(ngals)
   agelower = fltarr(ngals)
   fehupper = fltarr(ngals)
   fehlower = fltarr(ngals)
   for i=0,ngals-1 do begin
       chisqarr = make_chisqarr_new(feh(i),age(i),sn(i))
       chisq_cube[*,*,i] = chisqarr
       print,'finished gal no. ',i
   endfor
   objnamestr={objname:objname}
   if keyword_set(choi) then outname='chisqcube_choi_sci.fits' else $
       outname='chisqcube_sdss_sci.fits'
   mwrfits,chisq_cube,outname,/create,/silent
   mwrfits,grid_feh,outname,/silent
   mwrfits,grid_age,outname,/silent
   mwrfits,objnamestr,outname,/silent

end
