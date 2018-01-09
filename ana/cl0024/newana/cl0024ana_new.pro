pro linearfit,x,a,f,pder
       f=a[0]+a[1]*x
       pder = [[replicate(1.,n_elements(x))],[x]]
end


pro cl0024ana_new,redoprepcl=redoprepcl,redoprepsdss=redoprepsdss
;main program for analysis section of the paper. combine previous plot_mass*.pros
;need the matched catalog from gals_character.pro
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;GETTING SCI DATA
   if keyword_set(redoprepcl) then begin
      ;data retreiving
      sci = mrdfits('/scr2/nichal/workspace2/sps_fit/data/all_cl0024/sps_fit.fits.gz',1,/silent)
      cat = mrdfits('/scr2/nichal/keck/deimos/Cl0024MS0451/Cl0024master.v7.fits.gz',1,/silent)
      file_matchedcat = '/scr2/nichal/workspace2/ana/cl0024/cl0024_matchedcat.sav'
      ;oii equivalent widths
      restore,'/scr2/nichal/workspace2/ana/cl0024/cl0024_allsci_oii_ew.sav';oii_ew_arr,oii_ew_err_ar
      
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
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      ;categorize the fit
      wgoodmem = where(sci.goodfit+sci.good eq 2 and sci.zfit ge 0.37 and sci.zfit le 4.2 and sci.logmstar gt 7., cgals)
      wbadmem  = where(sci.goodfit+sci.good eq 1 and sci.zfit ge 0.37 and sci.zfit le 4.2 and sci.feh ne -999. and sci.logmstar gt 7., cbadgals)
      wallmem = [wgoodmem,wbadmem]
      
      goodmem = bytarr(nsci)
      goodmem(wgoodmem) = 1
      badmem = bytarr(nsci)
      badmem(wbadmem) = 1
      allmem = bytarr(nsci)
      allmem(wallmem) = 1
      
      wpassivemem = where(oii_ew_arr(wallmem)+oii_ew_err_arr(wallmem) gt -5,cpassivemem)
      wpassivemem = wallmem(wpassivemem)
      passivemem = bytarr(nsci)
      passivemem(wpassivemem) = 1
      
      wgoodpassivemem = where(passivemem eq 1 and goodmem eq 1,cgoodpassivemem)
      wbadpassivemem = where(passivemem eq 1 and badmem eq 1,cbadpassivemem)
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
      ; only take the interested quiescent galaxies
      sci = sci(wpassivemem) 
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    
      ;fix the uncertainties
      ;new style - this will also create a grid of prob deposit in a file chisqcube_cl0024_sci.fits
      get_additional_uncertainties_cl0024,sci.feh,sci.age,sci.sn,sci.objname
      
      ;create probability distribution from chisq cube
      get_prob_dist,filein='chisqcube_cl0024_sci.fits',fileout='cl0024_feh_age_probdist.fits'
 
     ;save these shit to its own fits file
      probsci = mrdfits('cl0024_feh_age_probdist.fits',1)
      struct_add_field,sci,'ageupper',sci.age+probsci.dageupper
      struct_add_field,sci,'agelower',sci.age+probsci.dagelower
      struct_add_field,sci,'fehupper',sci.feh+probsci.dfehupper
      struct_add_field,sci,'fehlower',sci.feh+probsci.dfehlower

      mwrfits,sci,'sci_cl0024_ana.fits',/create,/silent
   endif
      
   sci = mrdfits('sci_cl0024_ana.fits',1)
   probsci = mrdfits('cl0024_feh_age_probdist.fits',1)
   ageform = (galage(sci.zfit,1000.)/1.e9-sci.age)>0. ;age of universe when it was formed
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Averaging
   massbins = [8.9,9.3,9.6,10.,10.25,10.55,11.]
   nbins = n_elements(massbins)-1
   ave_mass    = fltarr(nbins)
   ave_mass_dev= fltarr(nbins)
   ave_feh     = fltarr(nbins) 
   ave_feh_dev = fltarr(nbins)
   bndry_mass = fltarr(nbins*2)
   for i=0,nbins-1 do begin
      if i eq nbins-1 then msel = where(sci.logmstar gt massbins[i],cmsel) else $
      msel = where(sci.logmstar gt massbins[i] and sci.logmstar lt massbins[i+1],cmsel)
      meanmc,sci(msel).feh,probsci(msel).dfeh,probsci(msel).probdfeh,fehmean,sigmam,sigmad,sigmas;,/plot
      ;meanerr,sci(msel).feh,(sci(msel).fehupper-sci(msel).fehlower)*0.5,fehmean,sigmam,sigmad,sigmas		
      ave_feh(i) = fehmean
      ave_feh_dev(i) = sigmas
      ave_mass(i) = mean(sci(msel).logmstar)
      ave_mass_dev(i) = stdev(sci(msel).logmstar)
      bndry_mass[i*2:i*2+1] = [massbins[i],massbins[i+1]]
   endfor
   hifeh = interpol(ave_feh+ave_feh_dev,ave_mass,bndry_mass)
   lofeh = interpol(ave_feh-ave_feh_dev,ave_mass,bndry_mass)
   wtoolow = where(lofeh lt -1.,cwtoolow)
   if cwtoolow ge 1 then lofeh(wtoolow) = -1.

   ;stop
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Get SDSS data
   if keyword_set(redoprepsdss) then begin
      sdss = mrdfits('/scr2/nichal/workspace2/sps_fit/data/gallazzi_allmass2/sps_fit.fits.gz',1,/silent)
      nofit = where(sdss.feh eq -999 ,cnofit)
      if cnofit gt 0 then remove,nofit,sdss
      restore,'/scr2/nichal/workspace2/SDSSana/gallazzi_allmass2_Ha_ew.sav' ;ha_ew_arr, ha_ew_err_arr
      sdssem = where(ha_ew_arr lt -2.,csdssem,complement=noem,ncomplement=cnoem)
      
      ;only take the non-star forming pop
      sdss = sdss(noem)
    
      restore,'/scr2/nichal/workspace2/SDSSana/gallazzi_allmass2/match_gallazzi_allmass.sav'
      catsdss = cat(noem)
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    
      ;fix the uncertainties
      ;new style - this will also create a grid of prob deposit in a file chisqcube_cl0024_sci.fits
      get_additional_uncertainties_sdss,sdss.feh,sdss.age,sdss.sn,sdss.objname

      ;create probability distribution from chisq cube
      get_prob_dist,filein='chisqcube_sdss_sci.fits',fileout='sdss_feh_age_probdist.fits'

     ;save these shit to its own fits file
      probsdss = mrdfits('sdss_feh_age_probdist.fits',1)
      struct_add_field,sdss,'ageupper',sdss.age+probsdss.dageupper
      struct_add_field,sdss,'agelower',sdss.age+probsdss.dagelower
      struct_add_field,sdss,'fehupper',sdss.feh+probsdss.dfehupper
      struct_add_field,sdss,'fehlower',sdss.feh+probsdss.dfehlower


      mwrfits,sdss,'sci_sdss_ana.fits',/create,/silent
      mwrfits,catsdss,'sci_sdss_ana.fits',/silent
   endif

   sdss = mrdfits('sci_sdss_ana.fits',1)
   catsdss = mrdfits('sci_sdss_ana.fits',2)
   probsdss = mrdfits('sdss_feh_age_probdist.fits',1)
   sdss_ageform = (galage(sdss.zfit,1000.)/1.e9-sdss.age)>0. ;age of universe when it was formed	
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Average the SDSS data
   massbins_sdss = [9,9.8,10.3,10.5,10.7,10.9,11.2,11.5]
   massbins_sdss = [9,9.4,9.8,10.1,10.4,10.7,10.9,11.2,11.5]
   nbins_sdss = n_elements(massbins_sdss)-1
   ave_mass_sdss    = fltarr(nbins_sdss)
   ave_mass_dev_sdss= fltarr(nbins_sdss)
   ave_feh_sdss     = fltarr(nbins_sdss)
   ave_feh_dev_sdss = fltarr(nbins_sdss)
   
   for i=0,nbins_sdss-1 do begin
      msel = where(sdss.logmstar gt massbins_sdss[i] and sdss.logmstar lt massbins_sdss[i+1],cmsel)
      meanmc,sdss(msel).feh,probsdss(msel).dfeh,probsdss(msel).probdfeh,fehmean,sigmam,sigmad,sigmas;,/plot
      ;meanerr,sdss(msel).feh,sdss(msel).feherr,fehmean,sigmam,sigmad,sigmas
      ave_feh_sdss(i) = fehmean
      ave_feh_dev_sdss(i) = sigmas
      ave_mass_sdss(i) = mean(sdss(msel).logmstar)
      ave_mass_dev_sdss(i) = stdev(sdss(msel).logmstar)
      ;if i eq nbins_sdss-1 then stop
   endfor
   bndry_mass_sdss = massbins_sdss
   hifeh_sdss = interpol(ave_feh_sdss+ave_feh_dev_sdss,ave_mass_sdss,bndry_mass_sdss)
   lofeh_sdss = interpol(ave_feh_sdss-ave_feh_dev_sdss,ave_mass_sdss,bndry_mass_sdss)
   mtoohi = where(bndry_mass_sdss gt 11.5,cmtoohi)
   if cmtoohi gt 0 then remove, mtoohi, bndry_mass_sdss, hifeh_sdss,lofeh_sdss
   
   ;Average the SDSS catalog data (measurements from G05)
   ave_mass_catsdss    = fltarr(nbins_sdss)
   ave_mass_dev_catsdss= fltarr(nbins_sdss)
   ave_feh_catsdss     = fltarr(nbins_sdss)
   ave_feh_dev_catsdss = fltarr(nbins_sdss)
   
   for i=0,nbins_sdss-1 do begin
   	msel = where(catsdss.logm50 gt massbins_sdss[i] and catsdss.logm50 lt massbins_sdss[i+1],cmsel)
   	meanerr,catsdss(msel).z50,(catsdss(msel).z84-catsdss(msel).z16)/2.,fehmean,sigmam,sigmad,sigmas
   	ave_feh_catsdss(i) = fehmean
   	ave_feh_dev_catsdss(i) = sigmas
   	ave_mass_catsdss(i) = mean(catsdss(msel).logm50)
   	ave_mass_dev_catsdss(i) = stdev(catsdss(msel).logm50)
   	;if i eq nbins_sdss-1 then stop
   endfor
   bndry_mass_catsdss = massbins_sdss
   hifeh_catsdss = interpol(ave_feh_catsdss+ave_feh_dev_catsdss,ave_mass_catsdss,bndry_mass_catsdss)
   lofeh_catsdss = interpol(ave_feh_catsdss-ave_feh_dev_catsdss,ave_mass_catsdss,bndry_mass_catsdss)
   mtoohi = where(bndry_mass_catsdss gt 11.5,cmtoohi)
   if cmtoohi gt 0 then remove, mtoohi, bndry_mass_catsdss, hifeh_catsdss,lofeh_catsdss
   
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;GETTING LITERATURE VALUES
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Get the average values measured in Gallazzi05
   g05_mass = [9.00,9.11,9.31,9.51,9.72,9.91,10.11,10.31,10.51,10.72,10.91,11.11,11.31,11.5];the first time is actually 8.91
   g05_feh  = [-0.6,-0.61,-0.65,-0.61,-0.52,-0.41,-0.23,-0.11,-0.01,0.04,0.07,0.10,0.12,0.13]
   ;below are the stars in the figure 8 of Gallazzi05
   g05_feherr = [0.62,0.56,0.59,0.55,0.47,0.43,0.35,0.31,0.27,0.25,0.22,0.21,0.2,0.2]/2.
   g05_fehlo = g05_feh-g05_feherr
   g05_fehhi = g05_feh+g05_feherr
   
   ;below are the diamonds in figure 8 of Gallazzi05
   ;g05_fehlo= [-1.11,-1.07,-1.1,-1.03,-0.97,-0.9,-0.8,-0.65,-0.41,-0.24,-0.14,-0.09,-0.06,-0.04]
   ;g05_fehhi= [0.0,0.0,-0.05,-0.01,0.05,0.09,0.14,0.17,0.20,0.22,0.24,0.25,0.26,0.28]
   toolow = where(g05_fehlo lt -1.,ctoolow)
   if ctoolow gt 0 then g05_fehlo(toolow) = -1

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Choi's Data
   z01 = {zlow:0.1,zhigh:0.2,mass:[9.9,10.2,10.4,10.7,11.0],Feh:[-0.05,-0.06,-0.01,-0.03,0.02],Feherr:[0.04,0.02,0.01,0.01,0.01]}
   z02 = {zlow:0.2,zhigh:0.3,mass:[10.2,10.5,10.7,11.0,11.3],Feh:[-0.08,-0.06,-0.03,-0.01,-0.05],Feherr:[0.04,0.02,0.01,0.01,0.02]}
   z03 = {zlow:0.3,zhigh:0.4,mass:[10.5,10.8,11.0,11.3],Feh:[-0.11,-0.05,-0.02,-0.03],Feherr:[0.03,0.01,0.01,0.02]}
   z04 = {zlow:0.4,zhigh:0.55,mass:[10.8,11.1,11.3],Feh:[-0.07,-0.04,-0.05],Feherr:[0.02,0.01,0.02]}
   z06 = {zlow:0.55,zhigh:0.7,mass:[10.9,11.0,11.3],Feh:[-0.15,-0.02,-0.05],Feherr:[0.07,0.03,0.03]}
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Xiangcheng's FIRE data
   readcol,'/scr2/nichal/workspace2/catalogs/xiangcheng_ma/mzr_z0pt8.txt',xma_mass08,xma_feh08
   readcol,'/scr2/nichal/workspace2/catalogs/xiangcheng_ma/mzr_z0.txt',xma_mass0,xma_feh0
   xma_feh08 = xma_feh08-0.2
   xma_feh0 = xma_feh0-0.2
   testma = 0
   if testma eq 1 then begin ; test for evolution
      good0 = where(xma_mass0 gt 9.3,cgood0)
      good08 = where(xma_mass08 gt 9.3,cgood08)
      ma_linfit_z0 = linfit(xma_mass0(good0)-10.,xma_feh0(good0),sigma=ma_linfit_z0_err)
      ma_linfit_z08 = linfit(xma_mass08(good08)-10.,xma_feh08(good08),sigma=ma_linfit_z08_err)
      ;slopes are not different
      ma_commonslope = wmean([ma_linfit_z0(1),ma_linfit_z08(1)],[ma_linfit_z0_err(1),ma_linfit_z08_err(1)])
      ;linfit with a fixed slope of 0.56684088 dex per log mass
      ma_linfit_fixslope_z0 = [ma_linfit_z0(0),ma_commonslope]
      yfitma0 = curvefit(xma_mass0(good0)-10.,xma_feh0(good0),fltarr(cgood0)+1,ma_linfit_fixslope_z0,ma_linfit_fixslope_z0_err,function_name='linearfit',fita=[1,0])
      ma_linfit_fixslope_z08 = [ma_linfit_z08(0),ma_commonslope]
      yfitma08 = curvefit(xma_mass08(good08)-10.,xma_feh08(good08),fltarr(cgood08)+1,ma_linfit_fixslope_z08,ma_linfit_fixslope_z08_err,function_name='linearfit',fita=[1,0])
      ma_zscore = (ma_linfit_fixslope_z0(0)-ma_linfit_fixslope_z08(0))/sqrt(ma_linfit_fixslope_z0_err(0)^2+ma_linfit_fixslope_z08_err(0)^2)
      print,'ma evolution = ',ma_linfit_fixslope_z0(0)-ma_linfit_fixslope_z08(0),'p/m',sqrt(ma_linfit_fixslope_z0_err(0)^2+ma_linfit_fixslope_z08_err(0)^2)
      print,'z score=',ma_zscore   
   ;  stop 
   endif
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;De Rossi 2017 (EAGLE)
   ;read off from Figure 5
   DeRossi_z0={mass:[9.15,9.48,9.81,10.15,10.5,10.72],feh:[-0.17,-0.08,0.05,0.12,0.17,0.31],feherr:[0.075,0.08,0.07,0.07,0.07,0.03]}
   DeRossi_z1={mass:[9.18,9.47,9.85,10.15,10.55],feh:[-0.39,-0.28,-0.09,0.08,0.17],feherr:[0.06,0.06,0.08,0.07,0.085]}
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Sybilska et al 2017 (hELENa, IFU from Sauron data)
   aa=read_csv('/scr2/nichal/workspace2/catalogs/sybilska.csv',n_table_header=1,header=header)
   Syb_z0 = {mass:aa.field1,feh:aa.field2}
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Lu 2014
   Lu_z0 = {mass:[8.21,8.61,9.0,9.41,9.8,10.2,10.61,11],feh:[-0.53,-0.43,-0.35,-0.26,-0.16,-0.098,-0.075,-0.09],$
            feherr:[0.13,0.12,0.11,0.12,0.11,0.10,0.095,0.095]}
   Lu_z1 = {mass:[8.21,8.61,9,9.4,9.8,10.2,10.6,11],feh:[-0.56,-0.44,-0.35,-0.25,-0.18,-0.13,-0.11,-0.11],$
             feherr:[0.14,0.13,0.11,0.11,0.11,0.10,0.10,0.09]}
   Lu_somerville_z0 = {mass:[8.25,8.65,9.05,9.45,9.85,10.25,10.64],feh:[-1.00,-0.82,-0.64,-0.46,-0.33,-0.15,-0.03],$
                      feherr:[0.06,0.1,0.1,0.1,0.1,0.10,0.095]}
   Lu_lu_z0 = {mass:[8.73,8.95,9.18,9.47,9.76,10.09,10.46,10.88],feh:[-.99,-0.84,-0.70,-0.52,-0.31,-0.06,0.18],$
              feherr:[0.1,0.1,0.1,0.1,0.1,0.1,0.10,0.095]}
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   ;Find the best fit MZR parameters
   ;1) functional form from Zahid13
   ;montecarlo
   zahid_mzr = fitzahidmc([sci.logmstar,sdss.logmstar],[sci.feh,sdss.feh],[[probsci.dfeh],[probsdss.dfeh]],$
               [[probsci.probdfeh],[probsdss.probdfeh]],zahid_mzr_sigma,chisq=chizahid)

   zahid_mzr_sdss = fitzahidmc(sdss.logmstar,sdss.feh,probsdss.dfeh,probsdss.probdfeh,zahid_mzr_sdss_sigma,chisq=chizahid_sdss) 

   zahid_mzr_cl = fitzahidmc(sci.logmstar,sci.feh,probsci.dfeh,probsci.probdfeh,zahid_mzr_cl_sigma,chisq=chizahid_cl)
 	
   ;2) 1 degree polynomial (linear equation - linfit)
   lin_mzr = linfitmc([sci.logmstar,sdss.logmstar]-10.,[sci.feh,sdss.feh],$
               [[probsci.dfeh],[probsdss.dfeh]],[[probsci.probdfeh],[probsdss.probdfeh]],lin_mzr_err,chisq=chilin)
 
   lin_mzr_sdss = linfitmc(sdss.logmstar-10.,sdss.feh,probsdss.dfeh,probsdss.probdfeh,lin_mzr_sdss_err,chisq=chilin_sdss,yfit=yfitsdss)

   lin_mzr_cl = linfitmc(sci.logmstar-10.,sci.feh,probsci.dfeh,probsci.probdfeh,lin_mzr_cl_err,chisq=chilin_cl,yfit=yfitcl)

   print,'best linear fit (combined):',sigfig([lin_mzr,lin_mzr_err],3)
   print,'best linear fit (sdss)    :',sigfig([lin_mzr_sdss,lin_mzr_sdss_err],3)
   print,'best linear fit (z=0.4)   :',sigfig([lin_mzr_cl,lin_mzr_cl_err],3)

   ;normal z test for slope and intercept
   ;slope
   zslope = (lin_mzr_sdss(1)-lin_mzr_cl(1))/sqrt(lin_mzr_sdss_err(1)^2+lin_mzr_cl_err(1)^2)
   print, 'normal z value for slope of SDSS and Cl0024 =',zslope
   ;The z value is ~-0.5 which is p=0.3. The slopes are not different. So can proceed to do intercept according to
   ;http://www.biostathandbook.com/ancova.html 
   ;find the weighted slope
   commonslope = (lin_mzr_cl(1)/lin_mzr_cl_err(1)^2+lin_mzr_sdss(1)/lin_mzr_sdss_err(1)^2)$
                 /(1./lin_mzr_cl_err(1)^2+1./lin_mzr_sdss_err(1)^2)
   commonslope_err = 1./sqrt(1./lin_mzr_cl_err(1)^2+1./lin_mzr_sdss_err(1)^2)
   commonslope_stdev = sqrt(((lin_mzr_cl(1)-commonslope)^2/lin_mzr_cl_err(1)^2+$
                             (lin_mzr_sdss(1)-commonslope)^2/lin_mzr_sdss_err(1)^2)/$
                            (1./lin_mzr_cl_err(1)^2+1./lin_mzr_sdss_err(1)^2)/2.)

   print, 'common slope',commonslope,commonslope_err,commonslope_stdev
   ;fit with common slope
   sdss_fixslope_pars = [lin_mzr_sdss(0),commonslope]
   yfit_sdss_fixslope = linfit_fixedslope_mc(sdss.logmstar-10.,sdss.feh,probsdss.dfeh,probsdss.probdfeh,$
                        sdss_fixslope_pars,sdss_fixslope_pars_err,chisq=chisq_sdss_fixslope)
   
   cl_fixslope_pars = [lin_mzr_cl(0),commonslope]
   yfit_cl_fixslope = linfit_fixedslope_mc(sci.logmstar-10.,sci.feh,probsci.dfeh,probsci.probdfeh,$
                        cl_fixslope_pars,cl_fixslope_pars_err,chisq=chisq_cl_fixslope)

   print, 'Fixed slope intercepts:'
   print, 'SDSS:',sdss_fixslope_pars(0),sdss_fixslope_pars_err(0)
   print, 'z=0.4:',cl_fixslope_pars(0),cl_fixslope_pars_err(0)
   ;mean of all mass
   meanmass=mean([sci.logmstar,sdss.logmstar]-10.)
   meanmass = 0. ;set to fit at 10^10 Msun
   print, 'mean mass:',meanmass
   sdss_feh_meanmass = commonslope*meanmass+sdss_fixslope_pars(0)
   sdss_feh_meanmass_err = sqrt(commonslope_err^2*meanmass^2+sdss_fixslope_pars_err(0)^2)
   cl_feh_meanmass = commonslope*meanmass+cl_fixslope_pars(0)
   cl_feh_meanmass_err = sqrt(commonslope_err^2*meanmass^2+cl_fixslope_pars_err(0)^2)
    
   ;z test comparing two means
   zscore = ((sdss_feh_meanmass-cl_feh_meanmass)-0.)/sqrt(sdss_feh_meanmass_err^2+cl_feh_meanmass_err^2)
   print, 'z score for intercepts:',zscore

   ;3) 2 degree polynomial
   poly_mzr =  polyfitmc([sci.logmstar,sdss.logmstar],[sci.feh,sdss.feh],2,[[probsci.dfeh],[probsdss.dfeh]],$
               [[probsci.probdfeh],[probsdss.probdfeh]],poly_mzr_err,chisq=chipoly)
   poly_mzr_sdss = polyfitmc(sdss.logmstar,sdss.feh,2,probsdss.dfeh,probsdss.probdfeh,poly_mzr_sdss_err,chisq=chipoly_sdss)
   poly_mzr_cl = polyfitmc(sci.logmstar,sci.feh,2,probsci.dfeh,probsci.probdfeh,poly_mzr_cl_err,chisq=chipoly_cl)

   ;4) linear fit with intrinsic scatter
   bestparam_sdss = intrinsic_scatter_linear_nongauss(sdss.logmstar-10.,sdss.feh,$
                    [sdss_fixslope_pars,0.06],probsdss.dfeh,probsdss.probdfeh,fita=[1,1,0])
   print, 'intercept, slope, int scatter(sdss)',sigfig(bestparam_sdss[0:2],3)
   bestparam_cl = intrinsic_scatter_linear_nongauss(sci.logmstar-10.,sci.feh,$
                  [cl_fixslope_pars,0.06],probsci.dfeh,probsci.probdfeh,fita=[1,1,0])
   print, 'intercept, slope, int scatter(z=0.4)',sigfig(bestparam_cl[0:2],3)
   ;stop
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Use Akaike's information criterion
   samplesize = n_elements([sci.logmstar,sdss.logmstar])
   k=2
   aic_lin = chilin+2.*k+(2.*k*(k+1.))/(samplesize-k-1) 
   k=3
   aic_poly= chipoly+2.*k+(2.*k*(k+1.))/(samplesize-k-1)
   aic_zahid = chizahid*(samplesize-3)+2.*k+(2.*k*(k+1.))/(samplesize-k-1)
   print, 'AIC linear = ', aic_lin
   print, 'AIC 2degre polynomial = ', aic_poly
   print, 'AIC Zahid =',aic_zahid
   
   print,'AIC Cl0024:linear(fixed slope),linear, polynomial'
   samplesize=n_Elements(sci.logmstar)
   k=1
   print,chisq_cl_fixslope+2.*k+(2.*k*(k+1.))/(samplesize-k-1)
   k=2
   print,chilin_cl+2.*k+(2.*k*(k+1.))/(samplesize-k-1)
   k=3
   print,chipoly_cl+2.*k+(2.*k*(k+1.))/(samplesize-k-1)
   
   print,'AIC SDSS:linear(fixed slope),linear, polynomial'
   samplesize=n_Elements(sdss)
   k=1
   print,chisq_sdss_fixslope+2.*k+(2.*k*(k+1.))/(samplesize-k-1)
   k=2
   print,chilin_sdss+2.*k+(2.*k*(k+1.))/(samplesize-k-1)
   k=3
   print,chipoly_sdss+2.*k+(2.*k*(k+1.))/(samplesize-k-1)
   
   ; ok, choose linear fit with fixed slope as the best fit according to AIC	
   mzr_mass = findgen(104)/40.+8.9
   mzr_zahid= zahid_mzr[0]-alog10(1.+10.^((mzr_mass-zahid_mzr[1])*(-1.*zahid_mzr[2])))-8.9
   mzr_poly = lin_mzr[0]+lin_mzr[1]*(mzr_mass-10.)
   mzr_poly_sdss = sdss_fixslope_pars[0]+sdss_fixslope_pars[1]*(mzr_mass-10.)
   mzr_poly_cl = cl_fixslope_pars[0]+cl_fixslope_pars[1]*(mzr_mass-10.)
   save, zahid_mzr,poly_mzr,poly_mzr_err,lin_mzr,lin_mzr_err,filename='mzr_obs_param.sav'

   ;random
   ;mean age of each population
    meanmc,sci.age,probsci.dage,probsci.probdage,sciagemean,sciagesigmam,sciagesigmad,sciagesigmas
    meanmc,sdss.age,probsdss.dage,probsdss.probdage,sdssagemean,sdssagesigmam,sdssagesigmad,sdssagesigmas 
    uniagez04 = galage(mean(sci.zfit),1000)/1.e9
    uniagez0  = galage(mean(sdss.zfit),1000)/1.e9
    sciageform = uniagez04-sciagemean
    sdssageform = uniagez0-sdssagemean
    print, (sciageform-sdssageform)*0.01
    ;stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;/////////////////////////////////////////////////////////////////////////////////////////////////////////////////;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;PLOTTING
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;misc
   set_plot,'ps'
   !p.multi = [0,1,1]
   !p.font = 0
   sunsym = sunsymbol()
   Delta = '!9'+string("104B)+'!x'
   alpha = '!9'+string("141B)+'!x'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   psname='SDSS_FeH_mass.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
   		xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      ;make the outline of plots
      xrange=[9,11.5]
      yrange=[-1.,0.3]
      plot,sdss.logmstar,sdss.feh,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
      
      ;shade the average region
      polyfill,[g05_mass,reverse(g05_mass)],[g05_fehhi,reverse(g05_fehlo)],color=fsc_color('pbg1')
      ;x=[bndry_mass,reverse(bndry_mass)]
      ;y=[hifeh,reverse(lofeh)]
      ;polyfill,x,y,color=fsc_color('org2')
      
      x=[bndry_mass_sdss,reverse(bndry_mass_sdss)]
      y=[hifeh_sdss,reverse(lofeh_sdss)]
      x(where(x eq 8.9)) = 9.
      polyfill,x,y,color=fsc_color('blu3'),/line_fill,orientation=45
      polyfill,x,y,color=fsc_color('blu3'),/line_fill,orientation=135
      rainbow_colors,n_colors=21
      ;draw axis
      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
      axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
      axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
      axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
      
      ;draw data points
      cgplot,sdss.logmstar,sdss.feh,psym=14,/overplot,color=40,symsize=1.3
      ;plot representative uncertainties
      lowermass = where(catsdss.logm50 lt 10.,complement=uppermass)
      cgerrplot,[-0.32,-0.32],[11.1+median(catsdss(lowermass).logm16-catsdss(lowermass).logm50),$
                               11.3+median(catsdss(uppermass).logm16-catsdss(uppermass).logm50)],$
                              [11.1+median(catsdss(lowermass).logm84-catsdss(lowermass).logm50),$
                               11.3+median(catsdss(uppermass).logm84-catsdss(uppermass).logm50)],$
                               /horizontal,color=40,thick=2
      cgerrplot,[11.1,11.3],[-0.32+median(sdss(lowermass).fehlower-sdss(lowermass).feh),$
                             -0.32+median(sdss(uppermass).fehlower-sdss(uppermass).feh)],$
                            [-0.32+median(sdss(lowermass).fehupper-sdss(lowermass).feh),$
                             -0.32+median(sdss(uppermass).fehupper-sdss(uppermass).feh)],$
                             color=40,thick=2

      ;Add DeRossi17 
      oploterror,derossi_z0.mass,derossi_z0.feh,derossi_z0.feherr,color=fsc_color('springgreen'),$
      linethick=2,errcolor=fsc_color('springgreen'),psym=1
      oplot,derossi_z0.mass,derossi_z0.feh,psym=cgsymcat(24),symsize=1.5,color=fsc_color('springgreen')
      
      ;Add Xiangcheng's data
      oplot,xma_mass0,xma_feh0,psym=cgsymcat(16),color=fsc_color('forestgreen'),symsize=1.2
      
       ;Add Choi's data
       oploterror,z01.mass,z01.feh,z01.feherr,color=fsc_color('red5'),linethick=2,errcolor=10
       oplot,z01.mass,z01.feh,psym=cgsymcat(46),color=fsc_color('red5'),symsize=2
      
      ;Add Sybilska2017
      oplot,syb_z0.mass,syb_z0.feh,color=fsc_color('maroon'),thick=3,linestyle=5
      
      ;Add best fitted line
      oplot,mzr_mass,mzr_poly_sdss,color=0,thick=3 ;linear function
      oplot,mzr_mass,zahid_mzr[0]-alog10(1.+10.^((mzr_mass-zahid_mzr[1])*(-1.*zahid_mzr[2])))-8.9,color=0,linestyle=3,thick=3
      oplot,mzr_mass,poly_mzr_sdss[0]+poly_mzr_sdss[1]*mzr_mass+poly_mzr_sdss[2]*mzr_mass^2,color=0,linestyle=3,thick=3
       ;Labelling
      cglegend,title=['Quiescent SDSS','measured in this work','Gallazzi et al. 2005','Choi et al. 2014','Sybilska et al. 2017'],psym=[14,0,15,46,0],location=[10.6,-0.5],box=0,charsize=0.8,/data,length=0,vspace=1.25,color=['navy','navy','pbg1','red5','maroon'],symsize=1.2

      cglegend,title=['Ma et al. 2016','De Rossi et al. 2017'],psym=[16,24],location=[10.6,-0.87],box=0,charsize=0.8,/data,length=0,vspace=1.25,color=['forestgreen','springgreen'],symsize=1.2

      oplot,[10.53,10.67],[-0.8,-0.8],linestyle=2,color=fsc_color('maroon'),thick=2
      oplot,[10.6],[-0.645],psym=cgsymcat(15),color=fsc_Color('pbg1'),symsize=1.3
      xyouts,10.1,-0.52,'observations:',charsize=0.8
      xyouts,10.1,-0.88,'simulations:',charsize=0.8
      xyouts,9.1,0.2,'z~0 galaxies',charsize=1.
      ;bracket,10.45,-0.8,10.5,-0.5,/left
      ;bracket,10.45,-0.95,10.5,-0.85,/left
   device,/close
   stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   psname='SDSS_FeH_mass_proposal.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
   	xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      ;make the outline of plots
      xrange=[9,11.5]
      yrange=[-1.,0.3]
      plot,sdss.logmstar,sdss.feh,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
      
      ;shade the average region
      polyfill,[g05_mass,reverse(g05_mass)],[g05_fehhi,reverse(g05_fehlo)],color=fsc_color('thistle')
      
      rainbow_colors,n_colors=21
      x=[bndry_mass_sdss,reverse(bndry_mass_sdss)]
      y=[hifeh_sdss,reverse(lofeh_sdss)]
      x(where(x eq 8.9)) = 9.
      polyfill,x,y,color=60,/line_fill,orientation=45
      polyfill,x,y,color=60,/line_fill,orientation=135
      
      ;draw axis
      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
      axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
      axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
      axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
      
      ;draw data points
      cgplot,sdss.logmstar,sdss.feh,psym=14,/overplot,color=40,symsize=1.3
      ;cgplot,sdss_em.logmstar,sdss_em.feh,psym=14,/overplot,color=40,symsize=1.3
      ;plot representative uncertainties
      lowermass = where(catsdss.logm50 lt 10.,complement=uppermass)
      cgerrplot,[-0.32,-0.32],[11.1+median(catsdss(lowermass).logm16-catsdss(lowermass).logm50),$
                               11.3+median(catsdss(uppermass).logm16-catsdss(uppermass).logm50)],$
                              [11.1+median(catsdss(lowermass).logm84-catsdss(lowermass).logm50),$
                               11.3+median(catsdss(uppermass).logm84-catsdss(uppermass).logm50)],$
                               /horizontal,color=40,thick=2
      cgerrplot,[11.1,11.3],[-0.32+median(sdss(lowermass).fehlower-sdss(lowermass).feh),$
                             -0.32+median(sdss(uppermass).fehlower-sdss(uppermass).feh)],$
                            [-0.32+median(sdss(lowermass).fehupper-sdss(lowermass).feh),$
                             -0.32+median(sdss(uppermass).fehupper-sdss(uppermass).feh)],$
                             color=40,thick=2
      ;oploterror,[11.1,11.3],[-0.32,-0.32],[median(0.5*catsdss(lowermass).logm84-0.5*catsdss(lowermass).logm16),$
      ;          median(0.5*catsdss(uppermass).logm84-0.5*catsdss(uppermass).logm16)],$
      ;          [0.5*median(sdss(lowermass).fehupper-sdss(lowermass).fehlower),$
      ;          0.5*median(sdss(uppermass).fehupper-sdss(uppermass).fehlower)],color=40,psym=3,errcolor=40,errthick=2
      
      ;Add Choi's data
      oploterror,z01.mass,z01.feh,z01.feherr,color=10,linethick=2,errcolor=10
      oplot,z01.mass,z01.feh,psym=cgsymcat(46),color=10,symsize=2
      
      ;Add Sybilska2017
      oplot,syb_z0.mass,syb_z0.feh,color=10,thick=3,linestyle=5
      
      ;Add best fitted line
      oplot,mzr_mass,mzr_poly_sdss,color=0,thick=3 ;linear function
      ;Labelling
      cglegend,title=['Leethochawalit et al. (in prep)','Gallazzi et al. 2005','Choi et al. 2014','Sybilska et al. 2017'],psym=[14,15,46,0],location=[10.3,-0.65],box=0,charsize=0.8,/data,length=0,vspace=1.25,color=[40,30,10,10],symsize=1.2
      oplot,[10.22,10.37],[-0.87,-0.87],linestyle=2,color=10,thick=2
      xyouts,9.1,0.2,'z~0 galaxies',charsize=1.
      ;bracket,10.45,-0.8,10.5,-0.5,/left
      ;bracket,10.45,-0.95,10.5,-0.85,/left
   device,/close
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
   cleanplot,/silent
   !p.font=0
   psname='Cl0024_FeH_mass.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      ;make the outline of plots
      xrange=[8.9,11.5]
      yrange=[-1.,0.3]
      plot,sci.logmstar,sci.feh,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
      ;shade the average region
      ;colorarr=[10,50,203,230,254]
      colorarr= fsc_color(['royalblue','darkorchid','deeppink','maroon','red8'])
      x=[bndry_mass,reverse(bndry_mass)]
      y=[hifeh,reverse(lofeh)]
      polyfill,x,y,color=fsc_color('rose')
      
      x=[bndry_mass_sdss,reverse(bndry_mass_sdss)]
      y=[hifeh_sdss,reverse(lofeh_sdss)]
      x(where(x eq 9)) = 8.9
      polyfill,x,y,color=fsc_color('blu3'),/line_fill,orientation=45
      polyfill,x,y,color=fsc_color('blu3'),/line_fill,orientation=135
      
      ;draw axis
      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
      axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
      axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
      axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
      
      ;Add Xiangcheng's data
      oplot,xma_mass0,xma_feh0,psym=cgsymcat(16),color=colorarr(0),symsize=1.2
      oplot,xma_mass08,xma_feh08,psym=cgsymcat(16),color=colorarr(4),symsize=1.2
      
      ;Add Munoz15
      
      ;Add DeRossi17 
      ;oploterror,derossi_z0.mass,derossi_z0.feh,derossi_z0.feherr,color=colorarr(0),linethick=2,errcolor=colorarr(0)
      ;oploterror,derossi_z1.mass,derossi_z1.feh,derossi_z1.feherr,color=colorarr(4),linethick=2,errcolor=colorarr(4)      
      ;oplot,derossi_z0.mass,derossi_z0.feh,psym=cgsymcat(24),color=colorarr(0),symsize=1.1
      ;oplot,derossi_z1.mass,derossi_z1.feh,psym=cgsymcat(24),color=colorarr(4),symsize=1.1
      
      ;draw data points
      cgerrplot,sci.logmstar,sci.fehlower,sci.fehupper,color='pink',thick=0.5
      oplot,sci.logmstar,sci.feh,psym=cgsymcat(14),color=colorarr(2),symsize=1.2
      oplot,sci.logmstar,sci.feh,psym=cgsymcat(4),color=fsc_color('maroon'),symsize=1.2
      
      ;Add Choi's data
      vsym,5,/fill,/star
      ;oploterror,z01.mass,z01.feh,z01.feherr,color=colorarr(1),linethick=2,errcolor=colorarr(1)
      oploterror,z03.mass,z03.feh,z03.feherr,color=fsc_color('org7'),linethick=2,errcolor=fsc_color('org7')
      ;oploterror,z06.mass,z06.feh,z06.feherr,color=colorarr(3),linethick=2,errcolor=colorarr(3)
      ;oplot,z01.mass,z01.feh,psym=cgsymcat(46),color=colorarr(1),symsize=2
      oplot,z03.mass,z03.feh,psym=cgsymcat(46),color=colorarr(2),symsize=1.7
      ;oplot,z06.mass,z06.feh,psym=cgsymcat(46),color=colorarr(3),symsize=2
      ;oplot,z01.mass,z01.feh,psym=cgsymcat(45),symsize=2
      oplot,z03.mass,z03.feh,psym=cgsymcat(45),symsize=1.7
      ;oplot,z06.mass,z06.feh,psym=cgsymcat(45),symsize=2
      
      ;Add Lu14
      ;oploterror,lu_z0.mass,lu_z0.feh,lu_z0.feherr,color=colorarr(0),linethick=2,errcolor=colorarr(0),errthick=2
      ;oploterror,lu_z1.mass,lu_z1.feh,lu_z1.feherr,color=colorarr(4),linethick=2,errcolor=colorarr(4)               
      oplot,lu_z0.mass,lu_z0.feh,psym=cgsymcat(24),color=colorarr(0),symsize=1.4
      ;oplot,lu_z1.mass,lu_z1.feh,psym=cgsymcat(24),color=colorarr(4),symsize=1.2
      
      
      ;add best fitted MZR curve
      ;oplot,mzr_mass,mzr_zahid,linestyle=2
     ; oplot,mzr_mass,mzr_poly,thick=2
      oplot,mzr_mass,mzr_poly_sdss,thick=3,color=fsc_color('navy'),linestyle=2
      oplot,mzr_mass,mzr_poly_cl,thick=3,color=fsc_color('org7'),linestyle=2
      ;Labelling
      ;al_legend,['z=0','z=[0.1,0.2]','z=[0.3,0.4]','z=[0.6,0.7]','z=[0.8,1]'],psym=15,color=colorarr,box=0,position=[10.6,-0.6]
      al_legend,['z=0','z=[0.3,0.4]','z=[0.8,1]'],psym=15,color=colorarr([0,2,4]),box=0,position=[10.6,-0.75]
      al_legend,['Current work','Choi et al. 2014','Ma et al. 2016','Lu et al. 2014'],color='black',psym=[14,46,16,24],position=[10.6,-0.47],box=0,charsize=0.9,symsize=[1.3,1.7,1.,1.2]		
      xyouts,9.,0.2,'z~0.4 galaxies',charsize=1.
   device,/close
stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;TEST KS TEST

   ;ks_test,sdss.logmstar,sdss.feh,sci.logmstar,sci.feh,D_ks,Neff_ks,pvalue   
   ; stop  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   psname='Cl0024_FeH_mass_proposal.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
            xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      ;make the outline of plots
      xrange=[8.9,11.5]
      yrange=[-1.,0.3]
      plot,sci.logmstar,sci.feh,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
      
      ;shade the average region
      ;colorarr=[10,50,203,230,254]
      colorarr= fsc_color(['darkorchid','royalblue','org5','maroon','red8'])
      x=[bndry_mass,reverse(bndry_mass)]
      y=[hifeh,reverse(lofeh)]
      polyfill,x,y,color=fsc_color('org2')
      
      x=[bndry_mass_sdss,reverse(bndry_mass_sdss)]
      y=[hifeh_sdss,reverse(lofeh_sdss)]
      x(where(x eq 9)) = 8.9
      polyfill,x,y,color=fsc_color('blu4'),/line_fill,orientation=45
      polyfill,x,y,color=fsc_color('blu4'),/line_fill,orientation=135
      
      ;draw axis
      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
      axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
      axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
      axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
      
      ;Add Munoz15
      
      ;Add DeRossi17 
      ;oploterror,derossi_z0.mass,derossi_z0.feh,derossi_z0.feherr,color=colorarr(0),linethick=2,errcolor=colorarr(0)
      ;oploterror,derossi_z1.mass,derossi_z1.feh,derossi_z1.feherr,color=colorarr(4),linethick=2,errcolor=colorarr(4)      
      ;oplot,derossi_z0.mass,derossi_z0.feh,psym=cgsymcat(24),color=colorarr(0),symsize=1.1
      ;oplot,derossi_z1.mass,derossi_z1.feh,psym=cgsymcat(24),color=colorarr(4),symsize=1.1
      ;draw data points
      cgerrplot,sci.logmstar,sci.fehlower,sci.fehupper,color='org3',thick=0.5
      oplot,sci.logmstar,sci.feh,psym=cgsymcat(14),color=colorarr(2),symsize=1.3
      oplot,sci.logmstar,sci.feh,psym=cgsymcat(4),color=fsc_color('org7'),symsize=1.3
      
      ;Add Xiangcheng's data
      ;oplot,xma_mass0,xma_feh0,psym=cgsymcat(16),color=colorarr(1),symsize=1.3
      ;oplot,xma_mass08,xma_feh08,psym=cgsymcat(16),color=colorarr(4),symsize=1.3
      ;Add Lu14
      ;oploterror,lu_z0.mass,lu_z0.feh,lu_z0.feherr,color=colorarr(0),linethick=2,errcolor=colorarr(0),errthick=2
      ;oplot,lu_z0.mass,lu_z0.feh,psym=cgsymcat(24),color=colorarr(0),symsize=1.2
      
      ;oploterror,lu_somerville_z0.mass,lu_somerville_z0.feh,lu_somerville_z0.feherr,color=colorarr(1),linethick=2,errcolor=colorarr(1),errthick=2
      ;oplot,lu_somerville_z0.mass,lu_somerville_z0.feh,psym=cgsymcat(24),color=colorarr(1),symsize=1.2
      ;oploterror,lu_lu_z0.mass,lu_lu_z0.feh,lu_lu_z0.feherr,color=colorarr(4),linethick=2,errcolor=colorarr(4),errthick=2
      ;oplot,lu_lu_z0.mass,lu_lu_z0.feh,psym=cgsymcat(24),color=colorarr(4),symsize=1.2
      
      ;Add Choi's data
      ;vsym,5,/fill,/star
      ;oploterror,z01.mass,z01.feh,z01.feherr,color=colorarr(1),linethick=2,errcolor=colorarr(1)
      oploterror,z03.mass,z03.feh,z03.feherr,color=fsc_color('org7'),linethick=2,errcolor=fsc_color('org7')
      ;oploterror,z06.mass,z06.feh,z06.feherr,color=colorarr(3),linethick=2,errcolor=colorarr(3)
      ;oplot,z01.mass,z01.feh,psym=cgsymcat(46),color=colorarr(1),symsize=2
      oplot,z03.mass,z03.feh,psym=cgsymcat(46),color=colorarr(2),symsize=2
      ;oplot,z06.mass,z06.feh,psym=cgsymcat(46),color=colorarr(3),symsize=2
      ;oplot,z01.mass,z01.feh,psym=cgsymcat(45),symsize=2
      oplot,z03.mass,z03.feh,psym=cgsymcat(45),symsize=2
      ;oplot,z06.mass,z06.feh,psym=cgsymcat(45),symsize=2
      
      ;add best fitted MZR curve
      ;oplot,mzr_mass,mzr_poly,thick=2
      oplot,mzr_mass,mzr_poly_sdss,thick=2,linestyle=2,color=fsc_color('navy')
      
      oplot,mzr_mass,mzr_poly_cl,thick=2,linestyle=2,color=fsc_color('red')
      ;oplot,mzr_mass,zahid_mzr_cl[0]-alog10(1.+10.^((mzr_mass-zahid_mzr_cl[1])*(-1.*zahid_mzr_cl[2])))-8.9,color=fsc_color('navy'),linestyle=1,thick=2
      ;oplot,mzr_mass,poly_mzr_cl[0]+poly_mzr_cl[1]*mzr_mass+poly_mzr_cl[2]*mzr_mass^2,color=fsc_color('maroon'),linestyle=2,thick=2
      ;oplot,[10.5,10.8],[-0.5,-0.5],thick=2
      ;oplot,[10.5,10.8],[-0.65,-0.65],color=fsc_color('navy'),linestyle=1,thick=2
      ;oplot,[10.5,10.8],[-0.8,-0.8],color=fsc_color('maroon'),linestyle=2,thick=2
      
      ;Labelling
      ;al_legend,['z~0 SDSS','z~0.4 Cl0024','combined'],psym=15,color=[colorarr[1],colorarr[2],fsc_color('black')],box=0,position=[10.6,-0.7]
      al_legend,['Leethochawalit et al. (in prep)','Choi et al. 2014'],color='black',psym=[14,46],position=[10.3,-0.6],box=0,charsize=0.9,symsize=[1.3,1.3]               
      al_legend,['z~0','z=[0.3,0.4]'],psym=15,color=colorarr([1,2]),box=0,position=[10.3,-0.75]
      ;al_legend,['z=0','z=[0.1,0.2]','z=[0.3,0.4]','z=[0.6,0.7]','z=[0.8,1]'],psym=15,color=colorarr,box=0,position=[10.6,-0.6]
      ;al_legend,['Current work','Lu+2014 (Croton)','Lu+2014 (Somerville)','Lu+2014 (Lu)'],color=[colorarr[2],colorarr[0],colorarr[1],colorarr[4]],psym=[14,24,24,24],position=[10.4,-0.6],box=0,charsize=0.9,symsize=[1.3,1.2,1.2,1.2]       
      
      xyouts,9.,0.2,'z~0.4 galaxies',charsize=1.
   device,/close
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Gas phase metal from Zahid2013 (12+log(O/H))
   gas_oh = [{redshift:0.08,z0:9.121,M0:8.999,gmma:0.85},{redshift:0.29,z0:9.130,M0:9.304,gmma:0.77},$
             {redshift:0.78,z0:9.161,M0:9.661,gmma:0.65},{redshift:1.4,z0:9.06,M0:9.6,gmma:0.7},$
             {redshift:2.26,z0:9.06,M0:9.7,gmma:0.6}]
   ;the equation is 12+log(O/H) = z0-log[1+(M*/M0)^-gmma]
   sun_oh = 8.8
;alpha_Fe = 0.12*mass-1.1 ;linear plot by eye to Choi14 Fig 8, mass is in log scale
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   psname='Formation_Redshift_cl0024.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
      	xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      zcat = [1000,2,1,0.7,0.4,0.]
      ;	zcat = [1000.,1.5,0.7,0]
      agecat = galage(zcat,1000.)/1.e9 ;age of universe at that redshif	
      plot,sci.logmstar,sci.feh,psym=1,xtitle='Log(M/M'+sunsym+')',ytitle='[Fe/H]',xrange=[9.,11.5],$
          xstyle=1,yrange=[-0.8,0.2],/nodata
      
      nage = n_Elements(agecat)-1
      rainbow_colors
      zcolor=reverse(fix((findgen(nage)+1.)/nage*254))
      ;loop over formation times
      for i=0,nage-1 do begin
      	selsdss = where(sdss_ageform gt agecat(i) and sdss_ageform le agecat(i+1),cselsdss)
      	if cselsdss gt 0 then begin
      		cgplot,sdss(selsdss).logmstar,sdss(selsdss).feh,psym=16,/overplot,color=zcolor(i),symsize=0.7
      	endif
      endfor
      for i=0,nage-1 do begin
              sel = where(ageform gt agecat(i) and ageform le agecat(i+1), csel)
              if csel gt 0 then begin
                cgerrplot,sci(sel).logmstar,sci(sel).fehlower,sci(sel).fehupper,color=zcolor(i)-10,thick=0.5
      		cgplot,sci(sel).logmstar,sci(sel).feh,psym=14,/overplot,color=zcolor(i),symsize=1.3
      		cgplot,sci(sel).logmstar,sci(sel).feh,psym=4,/overplot,color='darkgray',symsize=1.3
              endif
      endfor
      ;add gas phase MZR
      marr = findgen(101)/40.+9. ;logmass from 9 to 11.5
      a_fe = 0.12*marr-1.1
      for i=0,n_Elements(gas_oh)-1 do begin
      	o_h = gas_oh(i).z0-alog10(1.+10.^((marr-gas_oh(i).m0)*(-1.)*gas_oh(i).gmma))-sun_oh;-a_fe
      	colornowi = value_locate(zcat,gas_oh(i).redshift)
      	colornow = zcolor(colornowi)
      	oplot,marr,o_h,color=colornow,thick=0.5
      endfor
      oplot,mzr_mass,mzr_poly_sdss,thick=2,linestyle=2
      oplot,mzr_mass,mzr_poly_cl,thick=2,linestyle=2
      xyouts,[11.2,11.2],[0.,0.13],['z~0.4','z~0'],charsize=1
      ;Labelling
      zarr_str = strarr(n_elements(zcat)-1)
      for nz=0,n_elements(zcat)-2 do zarr_Str[nz]=strtrim(string(zcat[nz],format='(F3.1)'),2)+$
                                    '<z$\tex_{form}$<'+strtrim(string(zcat[nz+1],format='(F3.1)'),2)
      zarr_str(0) = 'z$\tex_{form}$>'+strtrim(string(zcat[1],format='(F3.1)'),2)
      al_Legend,zarr_str,psym=15,color=zcolor,box=0,thick=2,charsize=1,symsize=1.5,/right,/bottom,font=0
      al_Legend,['SDSS subsample','Cl0024 z~0.4'],psym=[16,14],symsize=[0.5,1.3],color=0,box=0,thick=2,$
                 charsize=1,position=[10.6,-0.3],font=0
      xyouts,9.1,0.1,'(a)',charsize=1.2
   device,/close
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   psname='FEH_deviation_obs.eps'
   fehideal = cl_fixslope_pars[0]+(sci.logmstar-10.)*cl_fixslope_pars[1]
   sdss_fehideal = sdss_fixslope_pars[0]+(sdss.logmstar-10.)*sdss_fixslope_pars[1]
   x=[sdss_ageform,ageform]
   y=[sdss.feh-sdss_fehideal,sci.feh-fehideal]
   dx = 0.5*[sdss.ageupper-sdss.agelower,sci.ageupper-sci.agelower]
   dy = 0.5*[sdss.fehupper-sdss.fehlower,sci.fehupper-sci.fehlower]
   dxarr = -1.*[[probsdss.dage],[probsci.dage]] ;since x axis is ageuni-agegal
   dyarr = [[probsdss.dfeh],[probsci.dfeh]]
   probdxdyarr = transpose([[[probsdss.probdfehdage]],[[probsci.probdfehdage]]],[1,0,2])
   feh_dev_par = linfitexymc(x,y,dx,dy,dxarr,dyarr,probdxdyarr,sigma_a_b)
   feh_dev_par_sdss = linfitexymc(sdss_ageform,sdss.feh-sdss_fehideal,0.5*[sdss.ageupper-sdss.agelower],0.5*[sdss.fehupper-sdss.fehlower],-1.*probsdss.dage,probsdss.dfeh,transpose(probsdss.probdfehdage,[1,0,2]),sigma_a_b_sdss)
   ;feh_dev_par_cl = linfitexymc(ageform,sci.feh-fehideal,0.5*[sci.ageupper-sci.agelower],0.5*[sci.fehupper-sci.fehlower],-1.*probsci.dage,probsci.dfeh,transpose(probsci.probdfehdage,[1,0,2]),sigma_a_b_cl,/plot)
   fitexy,ageform,sci.feh-fehideal,Acl,Bcl,x_sig=0.5*[sci.ageupper-sci.agelower],y_sig=0.5*[sci.fehupper-sci.fehlower],sigma_a_b_cl
   feh_dev_par_cl = [Acl,Bcl]
   print,'feh deviation'
   print,feh_dev_par,sigma_a_b
   print,'sdss:',feh_dev_par_sdss,sigma_a_b_sdss
   print,'cl:',feh_dev_par_cl,sigma_a_b_cl
   wslope = wmean([feh_dev_par_sdss[1],feh_dev_par_cl[1]],[sigma_a_b_sdss[1],sigma_a_b_cl[1]],error=wslope_err)
   
   set_plot,'ps'
   device, filename = psname,xsize = 15,ysize = 10, $
      xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      plot,sdss_ageform,sdss.feh-sdss_fehideal,psym=1,xtitle='Age of the universe at galaxy formation',ytitle=delta+'[Fe/H]',xrange=[0,14],xstyle=9,yrange=[-0.4,0.4],/nodata,position=[0.15,0.15,0.95,0.9]
      axis,xaxis=0,xstyle=1,xrange=[0,14]
      axis,xaxis=1,xticks=4,xtickv=[1.558,3.316, 5.903, 8.628,13.712],xtickn=['4','2','1','0.5','0'],xtitle='z'
      cgplot,sdss_ageform,sdss.feh-sdss_fehideal,psym=16,/overplot,symsize=0.5,color='ygb5'
      cgplot,ageform,sci.feh-fehideal,psym=14,/overplot,color='org4'
      oplot,[0,14],feh_dev_par(0)+feh_dev_par(1)*[0,14],thick=2
      ;oplot,[0,14],feh_dev_par_sdss(0)+feh_dev_par_sdss(1)*[0,14],thick=2,linestyle=2,color=fsc_color('blu5')
      ;oplot,[0,14],feh_dev_par_cl(0)+feh_dev_par_cl(1)*[0,14],thick=2,linestyle=2,color=fsc_color('org6')

      cgerrplot,[12.],-0.25+[median([sdss.fehlower-sdss.feh,sci.fehlower-sci.feh])],-0.25+[median([sdss.fehupper-sdss.feh,sci.fehupper-sci.feh])]
      cgerrplot,-0.25,12.+[median([sdss.agelower-sdss.age,sci.agelower-sci.age])],12.+[median([sdss.ageupper-sdss.age,sci.ageupper-sci.age])],/horizontal
      ;oploterror,[12.5],[-0.25],[median(dx)],[median(dy)]
      xyouts,0.5,0.32,'(c) real observations',charsize=1.2
      xyouts,1,-0.38,'older galaxies',charsize=0.8
      xyouts,10,-0.38,'younger galaxies',charsize=0.8
      arrow,0.9,-0.37,0.3,-0.37,/data
      arrow,13.1,-0.37,13.7,-0.37,/data
   device,/close
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;;Write table for latex
   make_catalog,sci.ra,sci.dec,sci.logmstar,sci.feh,sci.fehupper-sci.feh,sci.fehlower-sci.feh,$
                sci.age,sci.ageupper-sci.age,sci.agelower-sci.age,sci.snfit
   stop	
end

