pro match_oii_sci
;This program take the OII EW file whose entries match the all_cl0024/sps_fit.fits, find the match with the 
;final science sample in /ana/cl0024/newana/sci_cl0024_ana.fits and add a tag in the structure
;So hopefully running it once is enough
   sci = mrdfits('/scr2/nichal/workspace2/ana/cl0024/newana_after_referee_report/sci_cl0024_ana_afterrev.fits',1)
   sciall = mrdfits('/scr2/nichal/workspace2/sps_fit/data/all_cl0024/sps_fit.fits.gz',1,/silent)
   restore,'/scr2/nichal/workspace2/ana/cl0024/cl0024_allsci_oii_ew.sav';oii_ew_arr,oii_ew_err_ar
   ngals = n_Elements(sci)
   oii_sci = fltarr(ngals)
   oii_err_sci = fltarr(ngals)
   for i=0,ngals-1 do begin
      loc = where(sciall.objname eq sci[i].objname and sciall.mask eq sci[i].mask and $
                  sciall.slit eq sci[i].slit,cmatch)
      if cmatch ne 1 then stop
      oii_sci[i] = oii_ew_arr[loc]
      oii_err_sci[i] = oii_ew_err_arr[loc]
   endfor
   plothist,oii_sci[*,0],bin=1
   scitemp = sci[0]
   strout = create_struct(scitemp,'OIIew',-99.,'OIIew_err',-99.)
   strout = replicate(strout,ngals)
   struct_assign,sci,strout
   strout.oiiew = oii_sci
   strout.oiiew_err = oii_err_sci
   mwrfits,strout,'sci_cl0024_ana_afterrev.fits',/create,/silent

;Do SDSS
   sdss = mrdfits('/scr2/nichal/workspace2/ana/cl0024/newana/sci_sdss_ana.fits',1)
   sdssall = mrdfits('/scr2/nichal/workspace2/sps_fit/data/gallazzi_allmass2/sps_fit.fits.gz',1,/silent)
   restore,'/scr2/nichal/workspace2/SDSSana/gallazzi_allmass2_Ha_ew.sav' ;ha_ew_arr, ha_ew_err_arr
   sdssem = where(ha_ew_arr lt -2.,csdssem,complement=noem,ncomplement=cnoem)
   ngalsdss = n_Elements(sdss)
   if cnoem ne ngalsdss then stop,'something is wrong'
   ;only take the non-star forming pop
   ha_ew_arr = ha_ew_arr(noem)
   ha_ew_err_arr = ha_ew_err_arr(noem)
   strout = create_struct(sdss[0],'Haew',-99.,'Haew_err',-99.)
   strout = replicate(strout,ngalsdss)
   struct_assign,sdss,strout
   strout.haew = ha_ew_arr
   strout.haew_err = ha_ew_err_arr
   mwrfits,strout,'sci_sdss_ana_afterrev.fits',/create,/silent
   

end
