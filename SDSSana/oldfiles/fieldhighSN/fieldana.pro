pro fieldana
  ;READ DATA
  ;read science
  sci = mrdfits('/scr2/nichal/workspace2/sps_fit/data/field_highSN/sps_fit.fits.gz',1,/silent)
  ;read catalog
  cat = mrdfits('/scr2/nichal/workspace2/SDSSdata/field_highSN/Field_HighSN_nichaleet.fits',1,/silent)

  ;MATCHING
  ;match the science to catalog
  matcharr = lonarr(n_elements(sci))
  for i=0,n_elements(sci)-1 do begin
     match = where(sci[i].plate eq cat.plate and sci[i].mjd eq cat.mjd and sci[i].fiber eq cat.fiberid,cmatch)
     if cmatch ne 1 then stop
     matcharr[i] = match
  endfor
  ;rearrange the match
  for i=0, n_tags(cat)-1 do begin
     prematch = cat.(i)
     cat.(i)  = prematch(matcharr)
  endfor

  ;HISTOGRAMS/MAPS OF VARIOUS GENERAL INFO
  ;RA AND DEC
  set_plot,'ps'
  psname = 'RADEC.eps'
  device, filename = psname,xsize = 15,ysize = 12, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  plot, cat.ra, cat.dec,psym=4,xtitle='RA',ytitle='DEC'
  device,/close

  ;REDSHIFT,MASS DISTRIBUTION
  !p.multi = [0,2,1]
  psname = 'redshift_mass_hist.eps'
  device, filename = psname,xsize = 20,ysize = 8, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  plothist,cat.z,bin=0.01,xtitle='Redshift'
  plothist,cat.logmass,bin=0.1,xtitle='Log(Mass)',xrange=[9,13]
  device,/close

  ;MASS METALLICITY RELATION/ METALLICITY COMPARISON
  !p.multi = [0,3,2]
  !p.font  = 0 
  !p.charsize=2
  psname = 'field_MZR.eps'
  device, filename = psname,xsize = 30,ysize = 16, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  ;read gallazzi
  glz = mrdfits('/scr2/nichal/workspace2/gallazzi_data.fits',1,/silent)
  ;my measurement
  masserr = cat.logmass_err
  bad = where(masserr lt 0.)
  masserr(bad) = 0.
  ploterror,cat.logmass,sci.feh,masserr,sci.feherr,psym=1,xtitle='log M',ytitle='[Fe/H]',title='MZR (Early Form)',xrange=[9.5,12],xstyle=1,yrange=[-0.4,0.4],/nodata
  inmass = where(glz.logmass gt min(!x.crange))
  inz    = where(glz.zmin gt min(!y.crange))
  firstpoint = interpol(glz.z,glz.logmass,[min(!x.crange)])
  lastpoint  = interpol(glz.logmass,glz.zmin,[min(!y.crange)])
  cgcolorfill,[min(!x.crange),min(!x.crange),glz.logmass(inmass),reverse(glz.logmass(inz)),lastpoint],[min(!y.crange),firstpoint,glz.zmax(inmass),reverse(glz.zmin(inz)),min(!y.crange)],color=fsc_Color('cornsilk')
  oploterror,cat.logmass,sci.feh,masserr,sci.feherr,psym=1
  cgplot,cat.logmass,sci.feh,/overplot,psym=14
  oplot,glz.logmass,glz.z,color=fsc_color('chartreuse')
  oplot,glz.logmass,glz.zmin,psym=0,linestyle=2,color=fsc_color('chartreuse')
  oplot,glz.logmass,glz.zmax,psym=0,linestyle=2,color=fsc_color('chartreuse')
  ;find massbin everage
  ;eachbin has 20 member sort them by mass
  morder = sort(cat.logmass)
  for i=0,9 do begin
     in       = morder[i*20:(i+1)*20-1]
     lowmass = where(cat.logmass(in) lt 9.5,clowmass,complement=goodmass)
     if clowmass gt 0 then in = in(goodmass)
     massin   = cat.logmass(in)
     massinerr= cat.logmass_err(in)
     ifbad    = setintersection(bad,in)
     if ifbad[0] ne -1 then massinerr(where(massinerr lt 0)) = 1./0.
     massave  = wmean(massin,massinerr,/nan)
     fehin    = sci(in).feh
     fehinerr = sci(in).feherr
     badfeh   = where(fehinerr eq 0.,cbadfeh)
     if cbadfeh gt 0 then fehinerr(badfeh) = 1./0.
     fehave   = wmean(fehin,fehinerr,/nan)
     ;fehave   = mean(fehin)
     oploterror,[massave],[fehave],[stdev(massin)],[stdev(fehin)],errcolor=fsc_color('red')
     cgplot,[massave],[fehave],color=fsc_color('red'),/overplot,psym=14
     ;stop
  endfor

  ;SDSS measurement (early form)
  badmetal = where(cat.metallicity_err lt 0.)
  metalerr = cat.metallicity_err
  metalerr(badmetal) = 0.
  metal  = cat.metallicity ;zun is 0.019
  metal  = alog10(metal/0.019)
  ploterror,cat.logmass,metal,masserr,metalerr,psym=1,xtitle='log M',ytitle='SDSS ugriz metallicity (early form)',xrange=[9.5,12],xstyle=1
  cgplot,cat.logmass,metal,psym=14,/overplot
  cgplot,cat.logmass(bad),metal(bad),/overplot,psym=14,color=fsc_color('red')
  cgplot,cat.logmass(badmetal),metal(badmetal),/overplot,psym=14,color=fsc_color('purple')

  ;Compare Metallicity (early form)
  ;The differences
  metaldif = sci.feh-metal
  metaldiferr = sqrt(sci.feherr^2+metalerr^2)
  ploterror,cat.logmass,metaldif,masserr,metaldiferr,psym=1,xtitle='log M',ytitle='[Fe/H]-SDSS z(early form)',xrange=[9.5,12]
  cgplot,cat.logmass,metaldif,psym=14,/overplot
  cgplot,cat.logmass(badmetal),metaldif(badmetal),/overplot,psym=14,color=fsc_color('purple')

 ;my measurement 2
  masserr = cat.logmass_err_w
  bad = where(masserr lt 0.)
  masserr(bad) = 0.
  ploterror,cat.logmass_w,sci.feh,masserr,sci.feherr,psym=1,xtitle='log M',ytitle='[Fe/H]',title='field MZR (Late Form)',xrange=[9.5,12],xstyle=1,yrange=[-0.4,0.4]
  cgplot,cat.logmass_w,sci.feh,/overplot,psym=14
  cgplot,cat.logmass_w(bad),sci(bad).feh,/overplot,psym=14,color=fsc_color('red')

  ;SDSS measurement (late form)
  badmetal = where(cat.metallicity_err_w lt 0.)
  metalerr = cat.metallicity_err_w
  metalerr(badmetal) = 0.
  metal  = cat.metallicity_w ;zun is 0.019
  metal  = alog10(metal/0.019)
  ploterror,cat.logmass_w,metal,masserr,metalerr,psym=1,xtitle='log M',ytitle='SDSS ugriz metallicity(late form)',xrange=[9.5,12],xstyle=1
  cgplot,cat.logmass_w,metal,psym=14,/overplot
  cgplot,cat.logmass_w(bad),metal(bad),/overplot,psym=14,color=fsc_color('red')
  cgplot,cat.logmass_w(badmetal),metal(badmetal),/overplot,psym=14,color=fsc_color('purple')
  
  ;only 17 out of 200 have cat.logmass_w different from cat.logmass. All of the differences are less than 0.25 (in log space). The wide form have smaller mass than the early form.

  ;Compare Metallicity (wide form)
  ;The differences
  metaldif = sci.feh-metal
  metaldiferr = sqrt(sci.feherr^2+metalerr^2)
  ploterror,cat.logmass_w,metaldif,masserr,metaldiferr,psym=1,xtitle='log M',ytitle='[Fe/H]-SDSS z (late form)',xrange=[9.5,12]
  cgplot,cat.logmass_w,metaldif,psym=14,/overplot
  cgplot,cat.logmass_w(badmetal),metaldif(badmetal),/overplot,psym=14,color=fsc_color('purple')
  device,/close

  ;MZR BY CHI-SQ OF THE FIT
  !p.multi = [0,1,1]
  !p.charsize = 0
  psname = 'field_MZR_BY_CHISQ.eps'
  device, filename = psname,xsize = 10,ysize = 8, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color

  ;make the frame
  ploterror,cat.logmass,sci.feh,masserr,sci.feherr,psym=1,xtitle='log M',ytitle='[Fe/H]',title='MZR (Early Form)- Colored by Chisq',xrange=[9.5,12],xstyle=1,yrange=[-0.4,0.4],/nodata
  ;do the shading
  inmass = where(glz.logmass gt min(!x.crange))
  inz    = where(glz.zmin gt min(!y.crange))
  firstpoint = interpol(glz.z,glz.logmass,[min(!x.crange)])
  lastpoint  = interpol(glz.logmass,glz.zmin,[min(!y.crange)])
  cgcolorfill,[min(!x.crange),min(!x.crange),glz.logmass(inmass),reverse(glz.logmass(inz)),lastpoint],[min(!y.crange),firstpoint,glz.zmax(inmass),reverse(glz.zmin(inz)),min(!y.crange)],color=fsc_Color('corn silk')
  ;the real plots
  ;color bin by its chisq
  chisqbin = [0,5,10,20,40,100]
  color    = fsc_color(['navy','dodger blue','pink','deep pink','red'])
  for cb=0, n_Elements(chisqbin)-2 do begin
     inchi = where(sci.chisq gt chisqbin(cb) and sci.chisq lt chisqbin(cb+1),cinchi)
     if cinchi le 0 then stop,'PICK A NEW CHISQ BIN'
     oploterror,cat.logmass(inchi),sci(inchi).feh,masserr(inchi),sci(inchi).feherr,psym=1,errcolor=color(cb)
     cgplot,cat.logmass(inchi),sci(inchi).feh,/overplot,psym=14,color=color(cb)
     oplot,glz.logmass,glz.z,color=fsc_color('chartreuse')
     oplot,glz.logmass,glz.zmin,psym=0,linestyle=2,color=fsc_color('chartreuse')
     oplot,glz.logmass,glz.zmax,psym=0,linestyle=2,color=fsc_color('chartreuse')
  endfor

  ;find massbin everage - eachbin has 20 member sort them by mass
  ;only those with chisq lt 10.
  morder = sort(cat.logmass)
  inow = 0
  while inow lt n_elements(morder)-1 do begin
     nnow = 0
     min      = []
     minerr   = []
     fehin    = []
     fehinerr = []
     while nnow lt 20 and inow lt n_Elements(morder) do begin
        if cat.logmass(morder(inow)) gt 9.5 and sci(morder(inow)).chisq lt 10. then begin
           in = morder(inow)
           min = [min,cat.logmass(in)]
           minerr = [minerr,cat.logmass_err(in)]
           fehin  = [fehin,sci(in).feh]
           fehinerr = [fehinerr,sci(in).feherr]
           nnow += 1
        endif
        inow +=1
     endwhile
     badfeh = where(fehinerr eq 0.,cbadfeh)
     if cbadfeh gt 0 then fehinerr(badfeh) = 1./0.
     fehave   = wmean(fehin,fehinerr,/nan)
     massave  = wmean(min,minerr,/nan)
     if n_Elements(min) gt 1 then oploterror,[massave],[fehave],[stdev(min)],[stdev(fehin)],errcolor=fsc_color('red')
     cgplot,[massave],[fehave],color=fsc_color('black'),/overplot,psym=14
  endwhile
  al_legend,['0-5','5-10','10-20','20-40','40-100'],psym=14,color=color,box=0,margin=-0.5,/bottom
  device,/close

  ;VDISPERSION AND MASS
  ;to check how reliable is their mass measurement
  !p.multi = [0,1,1]
  psname = 'vdispersion_mass.eps'
  device, filename = psname,xsize = 10,ysize = 8, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  ploterror,cat.logmass,sci.vdisp,masserr,sci.vdisperr,psym=1,xtitle='log M', ytitle='vdisp (km/s)',xrange=[9.5,12],title='field mass check'
  oplot,cat.logmass(bad),sci(bad).vdisp,color=fsc_color('red'),psym=1
  device,/close
  
  ;AGE COMPARISONS
  !p.multi = [0,1,1]
  psname = 'age.eps'
  device, filename = psname,xsize = 10,ysize = 8, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  plot,cat.t_age,sci.age,psym=1,xtitle='SDSS age',ytitle='sps_fit Age (GYR)',/nodata
  cgErrplot,sci.age,cat.T_age_min,cat.T_age_max,/horizontal
  cgErrplot,cat.T_age,sci.age-sci.ageerr,sci.age+sci.ageerr
  cgplot,cat.t_age,sci.age,psym=14,/overplot
  device,/close
  stop
end
