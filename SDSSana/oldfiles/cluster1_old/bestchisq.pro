pro bestchisq
  common sps_spec, sps, spsz, spsage
  if (size(sps))[1] eq 0 then spsspec = sps_interp(0.0, 5.0)
  
;read science
  dir = '/scr2/nichal/workspace2/sps_fit/data/'
  scifits = dir+'cluster1/sps_fit.fits.gz'
  clustername = 'cluster1'
  scienceall = mrdfits(scifits,1,/silent)
;make output directory
  outputdir = '/scr2/nichal/workspace2/SDSSana/chisq/'+clustername
;select only members
  member = where(scienceall.zfit gt 0.08 and scienceall.zfit lt 0.12,cmember)
  scienceall = scienceall(member)
;sort them by mass
  sort_pos = sort(scienceall.logmstar)
  scienceall = scienceall(sort_pos)
;start doing the chisq check
  n_obj = n_elements(scienceall)
  chisqarr = fltarr(n_obj)
  badclump = bytarr(n_obj)
  for i=0,n_obj-1 do begin
     science      = scienceall[i]
     feh    = science.feh
     feherr = science.feherr
     zfit   = science.zfit
     age    = science.age
     ageerr = science.ageerr
     vdisp  = science.vdisp
     vdisperr = science.vdisperr
     lmass    = science.logmstar
   ;below is the same as sps_fit::fit
     restlambda = science.lambda / (1d + science.zspec)
     znow = science.zspec
     reallambda = science.lambda
     nlambda = n_elements(reallambda)
     if min(science.dlam) gt 0.2 and max(science.dlam,/nan) lt 10.0 then begin
        dlam_all = science.dlam
     endif else begin
        specresfile = dir+'cluster1/specres_poly.sav'
        if file_test(specresfile) then begin
           restore, specresfile
           dlam_all = poly(science.lambda/1000 - 7.8, specres_poly) / 2.35
        endif else dlam_all = replicate(3.9/2.35, nlambda)
     endelse

     won = where(science.fitmask eq 1 and finite(science.contdiv) and finite(science.contdivivar) and science.contdivivar gt 0, con)
    
     xmp  = reallambda[won]
     ymp  = science.contdiv[won]
     dymp = (science.contdivivar[won])^(-0.5)
     dlam = dlam_all[won]
     a = [feh,age,vdisp,zfit]
     curspec = get_sps_obs(xmp,a,dlam)
     chisq = total((curspec-ymp)^2/dymp^2)
     chisqarr[i] = chisq/float(n_elements(curspec)-4.) 
     badclump[i] = (feh lt .538*lmass-6.311 and lmass lt 10.8) or (feh lt -0.6 and lmass gt 10.8)
  endfor
goodchi = chisqarr(where(badclump eq 0))
badchi  = chisqarr(where(badclump))
end
