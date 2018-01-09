pro cluster1ana
  ;READ DATA
  ;read science
  sci = mrdfits('/scr2/nichal/workspace2/sps_fit/data/cluster1/sps_fit.fits.gz',1,/silent)
  ;read catalog
  cat = mrdfits('/scr2/nichal/workspace2/SDSSdata/cluster1/cluster1_nichaleet.fits',1,/silent)

  ;MATCHING
  ;match the science to catalog
  matcharr = lonarr(n_elements(sci))
  for i=0,n_elements(sci)-1 do begin
     match = where(sci[i].plate eq cat.plate and sci[i].mjd eq cat.mjd and sci[i].fiber eq cat.fiberid,cmatch)
     if cmatch ne 1 then stop
     matcharr[i] = match
  endfor

  ;remove bad measurement
  nofit = where(sci.feh eq -999 ,cnofit)
  if cnofit gt 0 then remove,nofit,sci,matcharr

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
  !p.multi = [0,3,1]
  !p.font  = 0 
  !p.charsize=2
  psname = 'cluster1_MZR.eps'
  device, filename = psname,xsize = 30,ysize = 8, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  ;read gallazzi
  glz = mrdfits('/scr2/nichal/workspace2/gallazzi_data.fits',1,/silent)
  ;my measurement
  masserr = cat.logmass_err
  bad = where(masserr lt 0.)
  masserr(bad) = 0.
  ploterror,cat.logmass,sci.feh,masserr,sci.feherr,psym=1,xtitle='log M (Early Form)',ytitle='[Fe/H]',title='MZR',xrange=[9.5,12],xstyle=1,yrange=[-1.5,0.4],/nodata
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
       plot,cat.logmass,sci.feh,psym=1,xtitle='log M (Early Form)',ytitle='[Fe/H]',title='MZR',xrange=[9.5,12],xstyle=1,yrange=[-1.5,0.4],/nodata,/noerase


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
  device,/close

  ;MZR BY CHI-SQ OF THE FIT
  !p.multi = [0,1,1]
  !p.charsize = 0
  psname = 'cluster1_MZR_BY_CHISQ.eps'
  device, filename = psname,xsize = 10,ysize = 8, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color

  ;make the frame
  ploterror,cat.logmass,sci.feh,masserr,sci.feherr,psym=1,xtitle='log M',ytitle='[Fe/H]',title='MZR (Early Form)- Colored by Chisq',xrange=[9.5,12],xstyle=1,yrange=[-1.4,0.4],/nodata
  ;do the shading
  inmass = where(glz.logmass gt min(!x.crange))
  inz    = where(glz.zmin gt min(!y.crange))
  firstpoint = interpol(glz.z,glz.logmass,[min(!x.crange)])
  lastpoint  = interpol(glz.logmass,glz.zmin,[min(!y.crange)])
  cgcolorfill,[min(!x.crange),min(!x.crange),glz.logmass(inmass),reverse(glz.logmass(inz)),lastpoint],[min(!y.crange),firstpoint,glz.zmax(inmass),reverse(glz.zmin(inz)),min(!y.crange)],color=fsc_Color('corn silk')
  ;the real plots
  ;color bin by its chisq
  chisqbin = [0,1,2,3,4,25]
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
  logmasstosort = cat.logmass
  if cnofit gt 0 then remove,findgen(cnofit)+n_Elements(sci),logmasstosort
  morder = sort(logmasstosort)
  inow = 0
  while inow lt n_elements(morder)-1-cnofit do begin
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
  al_legend,['0-1','1-2','2-3','3-4','4-25'],psym=14,color=color,box=0,margin=-0.5,/bottom
  device,/close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;MZR BY SIGNAL TO NOISE OF THE DATA
  !p.multi = [0,1,1]
  !p.charsize = 0
  psname = 'cluster1_MZR_BY_SN.eps'
  device, filename = psname,xsize = 10,ysize = 8, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color

  ;make the frame
  ploterror,cat.logmass,sci.feh,masserr,sci.feherr,psym=1,xtitle='log M',ytitle='[Fe/H]',title='MZR (Early Form)- Colored by SN',xrange=[9.5,12],xstyle=1,yrange=[-1.4,0.4],/nodata
  ;do the shading of galazzi
  inmass = where(glz.logmass gt min(!x.crange))
  inz    = where(glz.zmin gt min(!y.crange))
  firstpoint = interpol(glz.z,glz.logmass,[min(!x.crange)])
  lastpoint  = interpol(glz.logmass,glz.zmin,[min(!y.crange)])
  cgcolorfill,[min(!x.crange),min(!x.crange),glz.logmass(inmass),reverse(glz.logmass(inz)),lastpoint],[min(!y.crange),firstpoint,glz.zmax(inmass),reverse(glz.zmin(inz)),min(!y.crange)],color=fsc_Color('corn silk')
  ;the real plots
  ;color bin by its SN
  snbin = [0,10,15,20,25,40]
  color    = fsc_color(['navy','dodger blue','pink','deep pink','red'])
  for cb=0, n_Elements(snbin)-2 do begin
     insn = where(sci.sn gt snbin(cb) and sci.sn lt snbin(cb+1),cinsn)
     if cinsn le 0 then stop,'PICK A NEW SN BIN'
     oploterror,cat.logmass(insn),sci(insn).feh,masserr(insn),sci(insn).feherr,psym=1,errcolor=color(cb)
     cgplot,cat.logmass(insn),sci(insn).feh,/overplot,psym=14,color=color(cb)
     oplot,glz.logmass,glz.z,color=fsc_color('chartreuse')
     oplot,glz.logmass,glz.zmin,psym=0,linestyle=2,color=fsc_color('chartreuse')
     oplot,glz.logmass,glz.zmax,psym=0,linestyle=2,color=fsc_color('chartreuse')
  endfor

  ;find massbin everage - eachbin has 20 member sort them by mass
  ;only those with sn lt 10.
  logmasstosort = cat.logmass
  if cnofit gt 0 then remove,findgen(cnofit)+n_Elements(sci),logmasstosort
  morder = sort(logmasstosort)
  inow = 0
  while inow lt n_elements(morder)-1-cnofit do begin
     nnow = 0
     min      = []
     minerr   = []
     fehin    = []
     fehinerr = []
     while nnow lt 20 and inow lt n_Elements(morder) do begin
        if cat.logmass(morder(inow)) gt 9.5 and sci(morder(inow)).sn lt 10. then begin
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
  al_legend,['0-10','10-15','15-20','20-25','25-40'],psym=14,color=color,box=0,margin=-0.5,/bottom
  device,/close

  ;SN VS CHISQ
  !p.multi = [0,1,1]
  psname = 'sn_chisq.eps'
  device, filename = psname,xsize = 10,ysize = 8, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  plot,sci.sn,sci.chisq,psym=1,xtitle='SN',ytitle=cgsymbol('chi')+'!U2!N',/ylog,yrange=[0.8,25],ystyle=1
  device,/close
  

  ;VDISPERSION AND MASS
  ;to check how reliable is their mass measurement
  !p.multi = [0,1,1]
  psname = 'vdispersion_mass.eps'
  device, filename = psname,xsize = 10,ysize = 8, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  ploterror,cat.logmass,sci.vdisp,masserr,sci.vdisperr,psym=1,xtitle='log M', ytitle='vdisp (km/s)',xrange=[9.5,12],title='cluster1 mass check'
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
