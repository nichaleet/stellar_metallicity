pro maz_plot,masks=masks,all=all
;plot mass age metallicity measured from sps_fit

directory = '/scr2/nichal/workspace2/sps_fit/data/'
masks_path = file_search('/scr2/nichal/workspace2/sps_fit/data/','*',/test_directory)
nmasks = n_elements(masks_path)
masknames = replicate('mask',nmasks)

for i =0,nmasks-1 do masknames(i) = strmid(masks_path(i),1+strpos(masks_path(i),'/',/reverse_search),strlen(masks_path(i)))

if ~keyword_set(all) then begin
   masktemp=[]
   for ii=0,n_elements(masks)-1 do begin
      pos=where(strmatch(masknames,masks(ii),/fold_case),cpos)
      if cpos eq 0 then print,'Can''find mask ', masks(ii) else masktemp=[masktemp,masknames(pos)] 
   endfor
   masknames=masktemp
endif else masknames=masknames

for i=0,n_elements(masknames)-1 do begin
   sps_fit_all = mrdfits(directory+masknames(i)+'/sps_fit.fits.gz', 1, /silent)
   ;choose good spectrum
   goodspec= where(sps_fit.good eq 1)
   sps_fit = sps_fit_all(goodspec)

   ;group into z=0-0.5, z=0.5-1., z>1
   lowz = where(sps_fit.z lt 0.5 and sps_fit.z gt 0.)
   medz = where(sps_fit.z lt 1.  and sps_fit.z gt 0.5)
   highz = where(sps_fit.z gt 1. and sps_fit.z lt 2.)


   stop
endfor

end
