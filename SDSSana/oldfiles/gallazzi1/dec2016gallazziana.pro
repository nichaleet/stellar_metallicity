pro gallazziana
  ;READ DATA
  ;read science
  sci = mrdfits('/scr2/nichal/workspace2/sps_fit/data/gallazzi1/sps_fit.fits.gz',1,/silent)
  ;read fast measurement
  fast = []
  for i=0,19 do begin
     fastnow = readfastoutput('gallazzi1/gallazzi1_'+strtrim(string(fix(i*10)),2)+'_'+strtrim(string(fix(i*10+9)),2)+'.fout')
     fast = [fast,fastnow]
  endfor
  
  ;read catalog
  readcol,'/scr2/nichal/workspace2/SDSSdata/gallazzi1/z_mstar_all.dat',plate,mjd,fiber,logmasslow,logmass,logmasshigh,FeHlow,feh,fehhigh,logagelow,logage,logagehigh,comment='#'
  age = 10.^(logage-9.)
  agelow = 10.^(logagelow-9.)
  agehigh = 10.^(logagehigh-9.)
  catall = {plate:plate,mjd:mjd,fiber:fiber,logmasslow:logmasslow,logmass:logmass,logmasshigh:logmasshigh,fehlow:fehlow,feh:feh,fehhigh:fehhigh,agelow:agelow,age:age,agehigh:agehigh}

  ;remove bad measurement
  nofit = where(sci.feh eq -999 ,cnofit)
  if cnofit gt 0 then remove,nofit,sci,fast

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
  psname = 'redshift_mass_hist.eps'
  device, filename = psname,xsize = 20,ysize = 8, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  plothist,sci.zfit,bin=0.01,xtitle='Redshift'
  plothist,cat.logmass,bin=0.1,xtitle='Log(Mass)',xrange=[9,13]
  device,/close

  ;COMPARE THE MEASUREMENTS OF Z
  fehdiff = sci.feh-cat.feh
  fehdifflow = fehdiff-sqrt(sci.feherr^2+(cat.feh-cat.fehlow)^2)
  fehdiffhigh = fehdiff+sqrt(sci.feherr^2+(cat.feh-cat.fehhigh)^2)
  noem = where(sci.goodfit eq 1 and sci.feh ne 0.198, cnoem) ;no emission lines
  ferich = where(sci.goodfit eq 1 and sci.feh eq 0.198, cferich) ;no emission lines
  em = where(sci.goodfit eq 0, cem) ; with emission lines

  !p.multi=[0,3,2]
  !p.font=0
  !p.charsize=2
  psname = 'metal_comparison.eps'
  device,filename=psname,xsize=25,ysize=15,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  ;1)
  ploterror,cat.feh,sci.feh,sci.feherr,psym=1,xtitle='Gallazzi',ytitle='sps_fit',/nodata,yrange=[-1.5,0.5]
  oploterror,cat.feh(em),sci(em).feh,sci(em).feherr,psym=1,errcolor=fsc_color('black'),color=fsc_color('black')
  cgerrplot,sci(em).feh,cat.fehlow(em),cat.fehhigh(em),/horizontal,color='black'
  oploterror,cat.feh(ferich),sci(ferich).feh,sci(ferich).feherr,psym=1,errcolor=fsc_color('red'),color=fsc_color('red')
  cgerrplot,sci(ferich).feh,cat.fehlow(ferich),cat.fehhigh(ferich),/horizontal,color='red'
  oploterror,cat.feh(noem),sci(noem).feh,sci(noem).feherr,psym=1,errcolor=fsc_color('green'),color=fsc_color('green')
  cgerrplot,sci(noem).feh,cat.fehlow(noem),cat.fehhigh(noem),/horizontal,color='green'
  oplot,!x.crange,!x.crange,linestyle=2
  al_legend,['no EM lines','with EM lines'],psym=15,color=fsc_color(['green','black']),/bottom,/right,charsize=0
  ;2)
  ploterror,fast.feh,sci.feh,sci.feherr,psym=1,xtitle='FAST',ytitle='sps_fit',/nodata,yrange=[-1.5,0.5],xrange=[-1.5,0.5]
  oploterror,fast(em).feh,sci(em).feh,fast(em).feherr,sci(em).feherr,psym=1,errcolor=fsc_color('black'),color=fsc_color('black')
  oploterror,fast(ferich).feh,sci(ferich).feh,fast(ferich).feherr,sci(ferich).feherr,psym=1,errcolor=fsc_color('red'),color=fsc_color('red')
  oploterror,fast(noem).feh,sci(noem).feh,fast(noem).feherr,sci(noem).feherr,psym=1,errcolor=fsc_color('green'),color=fsc_color('green')
  oplot,!x.crange,!x.crange,linestyle=2
  ;3)
  plot,fast.feh,cat.feh,psym=1,xtitle='FAST',ytitle='Gallazzi',/nodata,yrange=[-1.5,0.5],xrange=[-1.5,0.5]
  cgerrplot,cat.feh(em),fast(em).feh-fast(em).feherr,fast(em).feh+fast(em).feherr,psym=1,color='black',/horizontal
  cgerrplot,fast(em).feh,cat.fehlow(em),cat.fehhigh(em),color='black'
  cgerrplot,cat.feh(ferich),fast(ferich).feh-fast(ferich).feherr,fast(ferich).feh+fast(ferich).feherr,psym=1,color='red',/horizontal
  cgerrplot,fast(ferich).feh,cat.fehlow(ferich),cat.fehhigh(ferich),color='red'
  cgerrplot,cat.feh(noem),fast(noem).feh-fast(noem).feherr,fast(noem).feh+fast(noem).feherr,psym=1,color='green',/horizontal
  cgerrplot,fast(noem).feh,cat.fehlow(noem),cat.fehhigh(noem),color='green'
  oplot,!x.crange,!x.crange,linestyle=2
  ;4)
  plot,cat.feh,fehdiff,xtitle='Gallazzi',ytitle='sps_fit-Gallazzi',/nodata,xrange=[-0.5,0.5];,yrange=[-1,0.5]
  cgerrplot,fehdiff(em),cat.fehlow(em),cat.fehhigh(em),/horizontal,color='black'
  cgerrplot,cat.feh(em),fehdifflow(em),fehdiffhigh(em),color='black'
  cgerrplot,fehdiff(ferich),cat.fehlow(ferich),cat.fehhigh(ferich),/horizontal,color='red'
  cgerrplot,cat.feh(ferich),fehdifflow(ferich),fehdiffhigh(ferich),color='red'
  cgerrplot,fehdiff(noem),cat.fehlow(noem),cat.fehhigh(noem),/horizontal,color='green'
  cgerrplot,cat.feh(noem),fehdifflow(noem),fehdiffhigh(noem),color='green'
  noemparams = linfit(cat.feh(noem),fehdiff(noem),measure_errors=0.5*(fehdiffhigh(noem)-fehdifflow(noem)),sigma=sigmanoem)
  oplot,!x.crange,noemparams[0]+noemparams[1]*!x.crange,color=fsc_color('darkgreen')
  ;5)
  fehdiff_fastsci = sci.feh-fast.feh
  fehdifferr_fastsci = sqrt(fast.feherr^2+sci.feherr^2)
  ploterror,fast.feh,fehdiff_fastsci,fast.feherr,fehdifferr_fastsci,xtitle='FAST',ytitle='sps_fit-FAST',/nodata
  oploterror,fast(em).feh,fehdiff_fastsci(em),fast(em).feherr,fehdifferr_fastsci(em),color=fsc_color('black'),errcolor=fsc_color('black'),psym=1
  oploterror,fast(ferich).feh,fehdiff_fastsci(ferich),fast(ferich).feherr,fehdifferr_fastsci(ferich),color=fsc_color('red'),errcolor=fsc_color('red'),psym=1
  oploterror,fast(noem).feh,fehdiff_fastsci(noem),fast(noem).feherr,fehdifferr_fastsci(noem),color=fsc_color('green'),errcolor=fsc_color('green'),psym=1
  noemparams = linfit(fast(noem).feh,fehdiff_fastsci(noem),measure_errors=fehdifferr_fastsci(noem),sigma=sigmanoem)
  oplot,!x.crange,noemparams[0]+noemparams[1]*!x.crange,color=fsc_color('darkgreen')

 ;6)
  fehdiff_fastgal = cat.feh-fast.feh
  fehdifflow_fastgal = fehdiff_fastgal-sqrt(fast.feherr^2+(cat.feh-cat.fehlow)^2)
  fehdiffhigh_fastgal = fehdiff_fastgal+sqrt(fast.feherr^2+(cat.feh-cat.fehhigh)^2)
  plot,fast.feh,fehdiff_fastgal,xtitle='FAST',ytitle='Gallazzi-FAST',/nodata
  cgerrplot,fehdiff_fastgal(em),fast(em).feh-fast(em).feherr,fast(em).feh+fast(em).feherr,/horizontal,color='black'
  cgerrplot,fast(em).feh,fehdifflow_fastgal(em),fehdiffhigh_fastgal(em),color='black'
  cgerrplot,fehdiff_fastgal(ferich),fast(ferich).feh-fast(ferich).feherr,fast(ferich).feh+fast(ferich).feherr,/horizontal,color='red'
  cgerrplot,fast(ferich).feh,fehdifflow_fastgal(ferich),fehdiffhigh_fastgal(ferich),color='red'
  cgerrplot,fehdiff_fastgal(noem),fast(noem).feh-fast(noem).feherr,fast(noem).feh+fast(noem).feherr,/horizontal,color='green'
  cgerrplot,fast(noem).feh,fehdifflow_fastgal(noem),fehdiffhigh_fastgal(noem),color='green'
  noemparams = linfit(fast(noem).feh,fehdiff_fastgal(noem),measure_errors=0.5*(fehdiffhigh_fastgal(noem)-fehdifflow_fastgal(noem)),sigma=sigmanoem)
  oplot,!x.crange,noemparams[0]+noemparams[1]*!x.crange,color=fsc_color('darkgreen')
  xyouts,0.5,0.97,'[Fe/H]',charsize=1.4,/normal
  
  ;plot,sci.feh,fehdiff,xtitle='sps_fit [Fe/H]',ytitle='sps_fit [Fe/H]-Gallazzi [Fe/H]',/nodata,xrange=[-0.5,0.5]
  ;cgerrplot,fehdiff(em),sci(em).feh-sci(em).feherr,sci(em).feh+sci(em).feherr,/horizontal,color='black'
  ;cgerrplot,sci(em).feh,fehdifflow(em),fehdiffhigh(em),color='black'
  ;cgerrplot,fehdiff(noem),sci(noem).feh-sci(noem).feherr,sci(noem).feh+sci(noem).feherr,/horizontal,color='green'
  ;cgerrplot,sci(noem).feh,fehdifflow(noem),fehdiffhigh(noem),color='green'
  device,/close

  !p.multi=[0,1,1]
  !p.font=0
  !p.charsize=0
  psname = 'metal_comparison_histogram.eps'
  device,filename=psname,xsize=10,ysize=10,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  plothist,cat.feh(noem),bin=0.05,ytitle='Galaxies with no EM',xtitle='[Fe/H]'
  plothist,sci(noem).feh,/overplot,color=fsc_color('blue'),bin=0.05
  plothist,cat.feh(noem),/overplot,color=fsc_color('red'),bin=0.05
  plothist,fast(noem).feh,/overplot,color=fsc_color('green'),bin=0.05
  al_legend,['sps_fit','Gallazzi','FAST'],psym=15,color=fsc_color(['blue','red','green'])
  device,/close





;COMPARE THE MEASUREMENTS OF AGE
  agediff = sci.age-cat.age
  agedifflow = agediff-sqrt(sci.ageerr^2+(cat.age-cat.agelow)^2)
  agediffhigh = agediff+sqrt(sci.ageerr^2+(cat.age-cat.agehigh)^2)
 
  !p.multi=[0,3,2]
  !p.font=0
  !p.charsize=2
  psname = 'age_comparison.eps'
  device,filename=psname,xsize=25,ysize=15,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  ;1)
  ploterror,cat.age,sci.age,sci.ageerr,psym=1,xtitle='Gallazzi',ytitle='sps_fit',/nodata
  oploterror,cat.age(em),sci(em).age,sci(em).ageerr,psym=1,errcolor=fsc_color('black'),color=fsc_color('black')
  cgerrplot,sci(em).age,cat.agelow(em),cat.agehigh(em),/horizontal,color='black'
  oploterror,cat.age(ferich),sci(ferich).age,sci(ferich).ageerr,psym=1,errcolor=fsc_color('red'),color=fsc_color('red')
  cgerrplot,sci(ferich).age,cat.agelow(ferich),cat.agehigh(ferich),/horizontal,color='red'
  oploterror,cat.age(noem),sci(noem).age,sci(noem).ageerr,psym=1,errcolor=fsc_color('green'),color=fsc_color('green')
  cgerrplot,sci(noem).age,cat.agelow(noem),cat.agehigh(noem),/horizontal,color='green'
  oplot,!x.crange,!x.crange,linestyle=2
  al_legend,['no EM lines','with EM lines'],psym=15,color=fsc_color(['green','black']),/top,/left,charsize=0
  ;2)
  ploterror,fast.age,sci.age,sci.ageerr,psym=1,xtitle='FAST',ytitle='sps_fit',/nodata
  oploterror,fast(em).age,sci(em).age,fast(em).ageerr,sci(em).ageerr,psym=1,errcolor=fsc_color('black'),color=fsc_color('black')
  oploterror,fast(ferich).age,sci(ferich).age,fast(ferich).ageerr,sci(ferich).ageerr,psym=1,errcolor=fsc_color('red'),color=fsc_color('red')
  oploterror,fast(noem).age,sci(noem).age,fast(noem).ageerr,sci(noem).ageerr,psym=1,errcolor=fsc_color('green'),color=fsc_color('green')
  oplot,!x.crange,!x.crange,linestyle=2
  ;3)
  plot,fast.age,cat.age,psym=1,xtitle='FAST',ytitle='Gallazzi',/nodata
  cgerrplot,cat.age(em),fast(em).age-fast(em).ageerr,fast(em).age+fast(em).ageerr,psym=1,color='black',/horizontal
  cgerrplot,fast(em).age,cat.agelow(em),cat.agehigh(em),color='black'
  cgerrplot,cat.age(ferich),fast(ferich).age-fast(ferich).ageerr,fast(ferich).age+fast(ferich).ageerr,psym=1,color='red',/horizontal
  cgerrplot,fast(ferich).age,cat.agelow(ferich),cat.agehigh(ferich),color='red'
  cgerrplot,cat.age(noem),fast(noem).age-fast(noem).ageerr,fast(noem).age+fast(noem).ageerr,psym=1,color='green',/horizontal
  cgerrplot,fast(noem).age,cat.agelow(noem),cat.agehigh(noem),color='green'
  oplot,!x.crange,!x.crange,linestyle=2
  ;4)
  plot,cat.age,agediff,xtitle='Gallazzi',ytitle='sps_fit-Gallazzi',/nodata
  cgerrplot,agediff(em),cat.agelow(em),cat.agehigh(em),/horizontal,color='black'
  cgerrplot,cat.age(em),agedifflow(em),agediffhigh(em),color='black'
  cgerrplot,agediff(ferich),cat.agelow(ferich),cat.agehigh(ferich),/horizontal,color='red'
  cgerrplot,cat.age(ferich),agedifflow(ferich),agediffhigh(ferich),color='red'
  cgerrplot,agediff(noem),cat.agelow(noem),cat.agehigh(noem),/horizontal,color='green'
  cgerrplot,cat.age(noem),agedifflow(noem),agediffhigh(noem),color='green'
  noemparams = linfit(cat.age(noem),agediff(noem),measure_errors=0.5*(agediffhigh(noem)-agedifflow(noem)),sigma=sigmanoem)
  oplot,!x.crange,noemparams[0]+noemparams[1]*!x.crange,color=fsc_color('darkgreen')
  ;5)
  agediff_fastsci = sci.age-fast.age
  agedifferr_fastsci = sqrt(fast.ageerr^2+sci.ageerr^2)
  ploterror,fast.age,agediff_fastsci,fast.ageerr,agedifferr_fastsci,xtitle='FAST',ytitle='sps_fit-FAST',/nodata
  oploterror,fast(em).age,agediff_fastsci(em),fast(em).ageerr,agedifferr_fastsci(em),color=fsc_color('black'),errcolor=fsc_color('black'),psym=1
  oploterror,fast(ferich).age,agediff_fastsci(ferich),fast(ferich).ageerr,agedifferr_fastsci(ferich),color=fsc_color('red'),errcolor=fsc_color('red'),psym=1
  oploterror,fast(noem).age,agediff_fastsci(noem),fast(noem).ageerr,agedifferr_fastsci(noem),color=fsc_color('green'),errcolor=fsc_color('green'),psym=1
  noemparams = linfit(fast(noem).age,agediff_fastsci(noem),measure_errors=agedifferr_fastsci(noem),sigma=sigmanoem)
  oplot,!x.crange,noemparams[0]+noemparams[1]*!x.crange,color=fsc_color('darkgreen')

 ;6)
  agediff_fastgal = cat.age-fast.age
  agedifflow_fastgal = agediff_fastgal-sqrt(fast.ageerr^2+(cat.age-cat.agelow)^2)
  agediffhigh_fastgal = agediff_fastgal+sqrt(fast.ageerr^2+(cat.age-cat.agehigh)^2)
  plot,fast.age,agediff_fastgal,xtitle='FAST',ytitle='Gallazzi-FAST',/nodata
  cgerrplot,agediff_fastgal(em),fast(em).age-fast(em).ageerr,fast(em).age+fast(em).ageerr,/horizontal,color='black'
  cgerrplot,fast(em).age,agedifflow_fastgal(em),agediffhigh_fastgal(em),color='black'
  cgerrplot,agediff_fastgal(ferich),fast(ferich).age-fast(ferich).ageerr,fast(ferich).age+fast(ferich).ageerr,/horizontal,color='red'
  cgerrplot,fast(ferich).age,agedifflow_fastgal(ferich),agediffhigh_fastgal(ferich),color='red'
  cgerrplot,agediff_fastgal(noem),fast(noem).age-fast(noem).ageerr,fast(noem).age+fast(noem).ageerr,/horizontal,color='green'
  cgerrplot,fast(noem).age,agedifflow_fastgal(noem),agediffhigh_fastgal(noem),color='green'
  noemparams = linfit(fast(noem).age,agediff_fastgal(noem),measure_errors=0.5*(agediffhigh_fastgal(noem)-agedifflow_fastgal(noem)),sigma=sigmanoem)
  oplot,!x.crange,noemparams[0]+noemparams[1]*!x.crange,color=fsc_color('darkgreen')
  xyouts,0.5,0.97,'AGE',charsize=1.4,/normal
  
  ;plot,sci.age,agediff,xtitle='sps_fit [Fe/H]',ytitle='sps_fit [Fe/H]-Gallazzi [Fe/H]',/nodata,xrange=[-0.5,0.5]
  ;cgerrplot,agediff(em),sci(em).age-sci(em).ageerr,sci(em).age+sci(em).ageerr,/horizontal,color='black'
  ;cgerrplot,sci(em).age,agedifflow(em),agediffhigh(em),color='black'
  ;cgerrplot,agediff(noem),sci(noem).age-sci(noem).ageerr,sci(noem).age+sci(noem).ageerr,/horizontal,color='green'
  ;cgerrplot,sci(noem).age,agedifflow(noem),agediffhigh(noem),color='green'
  device,/close

  !p.multi=[0,1,1]
  !p.font=0
  !p.charsize=0
  psname = 'age_comparison_histogram.eps'
  device,filename=psname,xsize=10,ysize=10,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  plothist,sci(noem).age,bin=1,ytitle='Galaxies with no EM',xtitle='[Fe/H]'
  plothist,sci(noem).age,/overplot,color=fsc_color('blue'),bin=1
  plothist,fast(noem).age,/overplot,color=fsc_color('green'),bin=1
  plothist,cat.age(noem),/overplot,color=fsc_color('red'),bin=1
  al_legend,['sps_fit','Gallazzi','FAST'],psym=15,color=fsc_color(['blue','red','green']),/top,/right
  device,/close

  !p.multi=[0,1,1]
  !p.font=0
  !p.charsize=0
  psname = 'age_metal_residuals.eps'
  device,filename=psname,xsize=10,ysize=10,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  plot,fehdiff,agediff,xtitle='Fe/H residuals',ytitle='AGE residuals',/nodata,xrange=[-0.8,0.8]
  cgerrplot,fehdiff(em),agedifflow(em),agediffhigh(em),color='black'
  cgerrplot,agediff(em),fehdifflow(em),fehdiffhigh(em),/horizontal,color='black'
  cgerrplot,fehdiff(ferich),agedifflow(ferich),agediffhigh(ferich),color='red'
  cgerrplot,agediff(ferich),fehdifflow(ferich),fehdiffhigh(ferich),/horizontal,color='red'
  cgerrplot,fehdiff(noem),agedifflow(noem),agediffhigh(noem),color='green'
  cgerrplot,agediff(noem),fehdifflow(noem),fehdiffhigh(noem),/horizontal,color='green'
  noemparams = linfit(fehdiff(noem),agediff(noem),measure_errors=0.5*(agediffhigh(noem)-agedifflow(noem)),sigma=sigmanoem)
  oplot,!x.crange,noemparams[0]+noemparams[1]*!x.crange,color=fsc_color('darkgreen')
  al_legend,['no EM lines','with EM lines'],psym=15,color=fsc_color(['green','black']),/top,/right,charsize=0
  device,/close
  stop
end
