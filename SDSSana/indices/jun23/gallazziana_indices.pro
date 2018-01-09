pro gallazziana_indices
  ;READ DATA
  ;read science
  sci = mrdfits('/scr2/nichal/workspace2/sps_fit/data/indices/gallazzi1/sps_fit.fits.gz',1,/silent)
  sciall = mrdfits('/scr2/nichal/workspace2/sps_fit/data/gallazzi1/sps_fit.fits.gz',1,/silent)
  
  ;read catalog
  readcol,'/scr2/nichal/workspace2/SDSSdata/gallazzi1/z_mstar_all.dat',plate,mjd,fiber,logmasslow,logmass,logmasshigh,FeHlow,feh,fehhigh,agelow,age,agehigh,comment='#'
  catall = {plate:plate,mjd:mjd,fiber:fiber,logmasslow:logmasslow,logmass:logmass,logmasshigh:logmasshigh,fehlow:fehlow,feh:feh,fehhigh:fehhigh,agelow:agelow,age:age,agehigh:agehigh}

  ;remove bad measurement
  nofit = where(sci.feh eq -999 ,cnofit)
  if cnofit gt 0 then remove,nofit,sci,matcharr

  ;MATCHING
  ;match the science to catalog
  matcharr = lonarr(n_elements(sci))
  for i=0,n_elements(sci)-1 do begin
     match = where(sci[i].plate eq catall.plate and sci[i].mjd eq catall.mjd and sci[i].fiber eq catall.fiber,cmatch)
     if cmatch ne 1 then stop
     matcharr[i] = match
  endfor

  ;get the match
  cat = {plate:plate(matcharr),mjd:mjd(matcharr),fiber:fiber(matcharr),logmasslow:logmasslow(matcharr),logmass:logmass(matcharr),logmasshigh:logmasshigh(matcharr),fehlow:fehlow(matcharr),feh:feh(matcharr),fehhigh:fehhigh(matcharr),agelow:agelow(matcharr),age:age(matcharr),agehigh:agehigh(matcharr)}


  ;REDSHIFT,MASS DISTRIBUTION
  set_plot,'ps'
  !p.multi = [0,2,1]
  psname = 'indices_redshift_mass_hist.eps'
  device, filename = psname,xsize = 20,ysize = 8, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  plothist,sci.zfit,bin=0.01,xtitle='Redshift'
  plothist,cat.logmass,bin=0.1,xtitle='Log(Mass)',xrange=[9,13]
  device,/close

  ;COMPARE THE MEASUREMENTS OF Z
  fehdiff = sci.feh-cat.feh
  fehdifflow = fehdiff-sqrt(sci.feherr^2+(cat.feh-cat.fehlow)^2)
  fehdiffhigh = fehdiff+sqrt(sci.feherr^2+(cat.feh-cat.fehhigh)^2)
  noem = where(sci.good eq 1 and sci.feh ne 0.198, cnoem) ;no emission lines
  ferich = where(sci.good eq 1 and sci.feh eq 0.198, cferich) ;no emission lines
  em = where(sci.good eq 0, cem) ; with emission lines

  !p.multi=[0,2,2]
  !p.font=0
  psname = 'indices_metal_comparison.eps'
  device,filename=psname,xsize=20,ysize=15,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  ploterror,cat.feh,sci.feh,sci.feherr,psym=1,xtitle='Gallazzi [Fe/H]',ytitle='sps_fit [Fe/H]',/nodata

  oploterror,cat.feh(em),sci(em).feh,sci(em).feherr,psym=1,errcolor=fsc_color('black'),color=fsc_color('black')
  cgerrplot,sci(em).feh,cat.fehlow(em),cat.fehhigh(em),/horizontal,color='black'

  oploterror,cat.feh(ferich),sci(ferich).feh,sci(ferich).feherr,psym=1,errcolor=fsc_color('red'),color=fsc_color('red')
  cgerrplot,sci(ferich).feh,cat.fehlow(ferich),cat.fehhigh(ferich),/horizontal,color='red'
  oploterror,cat.feh(noem),sci(noem).feh,sci(noem).feherr,psym=1,errcolor=fsc_color('green'),color=fsc_color('green')
  cgerrplot,sci(noem).feh,cat.fehlow(noem),cat.fehhigh(noem),/horizontal,color='green'
  oplot,!x.crange,!x.crange,linestyle=2
  al_legend,['no EM lines','with EM lines'],psym=15,color=fsc_color(['green','black']),/bottom,/right

  plot,cat.feh,fehdiff,xtitle='Gallazzi [Fe/H]',ytitle='sps_fit [Fe/H]-Gallazzi [Fe/H]',/nodata,xrange=[-0.5,0.5],yrange=[-1,0.5]
  cgerrplot,fehdiff(em),cat.fehlow(em),cat.fehhigh(em),/horizontal,color='black'
  cgerrplot,cat.feh(em),fehdifflow(em),fehdiffhigh(em),color='black'
  cgerrplot,fehdiff(noem),cat.fehlow(noem),cat.fehhigh(noem),/horizontal,color='green'
  cgerrplot,cat.feh(noem),fehdifflow(noem),fehdiffhigh(noem),color='green'
  noemparams = linfit(cat.feh(noem),fehdiff(noem),measure_errors=0.5*(fehdiffhigh(noem)-fehdifflow(noem)),sigma=sigmanoem)
  oplot,!x.crange,noemparams[0]+noemparams[1]*!x.crange,color=fsc_color('darkgreen')

  plot,sci.feh,fehdiff,xtitle='sps_fit [Fe/H]',ytitle='sps_fit [Fe/H]-Gallazzi [Fe/H]',/nodata,xrange=[-0.5,0.5]
  cgerrplot,fehdiff(em),sci(em).feh-sci(em).feherr,sci(em).feh+sci(em).feherr,/horizontal,color='black'
  cgerrplot,sci(em).feh,fehdifflow(em),fehdiffhigh(em),color='black'
  cgerrplot,fehdiff(noem),sci(noem).feh-sci(noem).feherr,sci(noem).feh+sci(noem).feherr,/horizontal,color='green'
  cgerrplot,sci(noem).feh,fehdifflow(noem),fehdiffhigh(noem),color='green'

  plothist,cat.feh(noem),bin=0.05,ytitle='Galaxies with no EM',xtitle='[Fe/H]'
  plothist,sci(noem).feh,/overplot,color=fsc_color('blue'),bin=0.05
  plothist,cat.feh(noem),/overplot,color=fsc_color('red'),bin=0.05
  al_legend,['sps_fit','Gallazzi'],psym=15,color=fsc_color(['blue','red'])
  device,/close

;COMPARE THE MEASUREMENTS OF AGE
  agediff = sci.age-cat.age
  agedifflow = agediff-sqrt(sci.ageerr^2+(cat.age-cat.agelow)^2)
  agediffhigh = agediff+sqrt(sci.ageerr^2+(cat.age-cat.agehigh)^2)
 
  !p.multi=[0,2,2]
  !p.font=0
  psname = 'indices_age_comparison.eps'
  device,filename=psname,xsize=20,ysize=15,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  ploterror,cat.age,sci.age,sci.ageerr,psym=1,xtitle='Gallazzi age',ytitle='sps_fit age',/nodata

  oploterror,cat.age(em),sci(em).age,sci(em).ageerr,psym=1,errcolor=fsc_color('black'),color=fsc_color('black')
  cgerrplot,sci(em).age,cat.agelow(em),cat.agehigh(em),/horizontal,color='black'

  oploterror,cat.age(ferich),sci(ferich).age,sci(ferich).ageerr,psym=1,errcolor=fsc_color('red'),color=fsc_color('red')
  cgerrplot,sci(ferich).age,cat.agelow(ferich),cat.agehigh(ferich),/horizontal,color='red'
  oploterror,cat.age(noem),sci(noem).age,sci(noem).ageerr,psym=1,errcolor=fsc_color('green'),color=fsc_color('green')
  cgerrplot,sci(noem).age,cat.agelow(noem),cat.agehigh(noem),/horizontal,color='green'
  oplot,!x.crange,!x.crange,linestyle=2
  al_legend,['no EM lines','with EM lines'],psym=15,color=fsc_color(['green','black'])

  plot,cat.age,agediff,xtitle='Gallazzi age',ytitle='sps_fit age-Gallazzi age',/nodata
  cgerrplot,agediff(em),cat.agelow(em),cat.agehigh(em),/horizontal,color='black'
  cgerrplot,cat.age(em),agedifflow(em),agediffhigh(em),color='black'
  cgerrplot,agediff(noem),cat.agelow(noem),cat.agehigh(noem),/horizontal,color='green'
  cgerrplot,cat.age(noem),agedifflow(noem),agediffhigh(noem),color='green'
  noemparams = linfit(cat.age(noem),agediff(noem),measure_errors=0.5*(agediffhigh(noem)-agedifflow(noem)),sigma=sigmanoem)
  oplot,!x.crange,noemparams[0]+noemparams[1]*!x.crange,color=fsc_color('darkgreen')

  plot,sci.age,agediff,xtitle='sps_fit age',ytitle='sps_fit age-Gallazzi age',/nodata
  cgerrplot,agediff(em),sci(em).age-sci(em).ageerr,sci(em).age+sci(em).ageerr,/horizontal,color='black'
  cgerrplot,sci(em).age,agedifflow(em),agediffhigh(em),color='black'
  cgerrplot,agediff(noem),sci(noem).age-sci(noem).ageerr,sci(noem).age+sci(noem).ageerr,/horizontal,color='green'
  cgerrplot,sci(noem).age,agedifflow(noem),agediffhigh(noem),color='green'

  plothist,sci(noem).age,bin=1,ytitle='Galaxies with no EM',xtitle='age'
  plothist,sci(noem).age,/overplot,color=fsc_color('blue'),bin=1
  plothist,cat.age(noem),/overplot,color=fsc_color('red'),bin=1
  al_legend,['sps_fit','Gallazzi'],psym=15,color=fsc_color(['blue','red'])
  device,/close

;COMPARE THE MEASUREMENTS OF Z between the two SPS_FIT
  fehdiff = sci.feh-sciall.feh
  fehdifflow = fehdiff-sqrt(sci.feherr^2+sciall.feherr^2)
  fehdiffhigh = fehdiff+sqrt(sci.feherr^2+sciall.feherr^2)
  noem = where(sci.good eq 1 and sci.feh ne 0.198 and sciall.feherr ne 0.198, cnoem) ;no emission lines
  ferich = where(sci.good eq 1 and sci.feh eq 0.198 or sciall.feh eq 0.198, cferich) ;no emission lines
  em = where(sci.good eq 0, cem) ; with emission lines

  !p.multi=[0,2,2]
  !p.font=0
  psname = 'sps_metal_comparison.eps'
  device,filename=psname,xsize=20,ysize=15,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  ploterror,sciall.feh,sci.feh,sci.feherr,psym=1,xtitle='all spec [Fe/H]',ytitle='indices region [Fe/H]',/nodata

  oploterror,sciall(em).feh,sci(em).feh,sciall(em).feherr,sci(em).feherr,psym=1,errcolor=fsc_color('black'),color=fsc_color('black')

  oploterror,sciall(ferich).feh,sci(ferich).feh,sciall(ferich).feherr,sci(ferich).feherr,psym=1,errcolor=fsc_color('red'),color=fsc_color('red')

  oploterror,sciall(noem).feh,sci(noem).feh,sciall(noem).feherr,sci(noem).feherr,psym=1,errcolor=fsc_color('green'),color=fsc_color('green')

  oplot,!x.crange,!x.crange,linestyle=2
  al_legend,['no EM lines','with EM lines'],psym=15,color=fsc_color(['green','black']),/bottom,/right

  plot,sciall.feh,fehdiff,xtitle='all spec [Fe/H]',ytitle='indices [Fe/H]-allspec [Fe/H]',/nodata,xrange=[-0.5,0.5]
  oploterror,sciall(em).feh,fehdiff(em),sciall(em).feherr,0.5*(fehdiffhigh(em)-fehdifflow(em)),color=fsc_color('black'),psym=1
  oploterror,sciall(noem).feh,fehdiff(noem),sciall(noem).feherr,0.5*(fehdiffhigh(noem)-fehdifflow(noem)),errcolor=fsc_color('green'),psym=1
  noemparams = linfit(sciall(noem).feh,fehdiff(noem),measure_errors=0.5*(fehdiffhigh(noem)-fehdifflow(noem)),sigma=sigmanoem)
  oplot,!x.crange,noemparams[0]+noemparams[1]*!x.crange,color=fsc_color('darkgreen')

  plot,sci.feh,fehdiff,xtitle='indices [Fe/H]',ytitle='indices [Fe/H]-all spec [Fe/H]',/nodata,xrange=[-0.5,0.5]
  oploterror,sci(em).feh,fehdiff(em),sci(em).feherr,0.5*(fehdiffhigh(em)-fehdifflow(em)),color=fsc_color('black'),psym=1
  oploterror,sci(noem).feh,fehdiff(noem),sci(noem).feherr,0.5*(fehdiffhigh(noem)-fehdifflow(noem)),errcolor=fsc_color('green'),psym=1


  plothist,sciall(noem).feh,bin=0.05,ytitle='Galaxies with no EM',xtitle='[Fe/H]'
  plothist,sci(noem).feh,/overplot,color=fsc_color('blue'),bin=0.05
  plothist,sciall(noem).feh,/overplot,color=fsc_color('red'),bin=0.05
  al_legend,['indices region','all spec'],psym=15,color=fsc_color(['blue','red'])
  device,/close

  stop
end
