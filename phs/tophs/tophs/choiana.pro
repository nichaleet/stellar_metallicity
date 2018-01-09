pro choiana
  ;READ DATA
  ;read science
  all = mrdfits('sps_fit_all.fits.gz',1,/silent)
  noMgb = mrdfits('sps_fit_noMgb.fits.gz',1,/silent)
  
  polycont = mrdfits('sps_fit_polycontn7.fits.gz',1,/silent)

  polycontv = mrdfits('sps_fit_polycontn4to7.fits.gz',1,/silent)

  ;COMPARE THE MEASUREMENTS OF Z
  fehdiff = all.feh-all.fehchoi
  fehdifferr = sqrt(all.feherr^2+all.fehchoierr^2)
  set_plot,'ps'
  !p.multi=[0,2,2]
  !p.font=0
  psname = 'choi_metal_comparison_all.eps'
  device,filename=psname,xsize=20,ysize=15,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  ploterror,all.fehchoi,all.feh,all.fehchoierr,all.feherr,psym=1,xtitle='Choi',ytitle='sps_fit'
  oplot,!x.crange,!x.crange,linestyle=2

  ploterror,all.fehchoi,fehdiff,all.fehchoierr,fehdifferr,xtitle='Choi',ytitle='sps_fit-Choi ',psym=1
  params = linfit(all.fehchoi,fehdiff,measure_errors=fehdifferr,sigma=sigmaparams)
  oplot,!x.crange,[-0.1,-0.1],linestyle=1,color=fsc_color('cyan')
  oplot,!x.crange,[0.1,0.1],linestyle=1,color=fsc_color('cyan')

  ploterror,all.MgFechoi,fehdiff,all.MgFechoierr,fehdifferr,xtitle='Choi [Mg/Fe]',ytitle='sps_fit-Choi ',psym=1

  plothist,fehdiff,bin=0.05,ytitle='Galaxies',xtitle='Choi-sps_fit'
  xyouts,0.5,0.97,'[Fe/H], All wavelengths',charsize=1.3,alignment=0.5,/normal
  device,/close

;COMPARE THE MEASUREMENTS OF AGE
  agediff = all.age-all.agechoi
  agedifferr = sqrt(all.ageerr^2+all.agechoierr^2)

  !p.multi=[0,2,2]
  !p.font=0
  psname = 'choi_age_comparison_all.eps'
  device,filename=psname,xsize=20,ysize=15,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  ploterror,all.agechoi,all.age,all.agechoierr,all.ageerr,psym=1,xtitle='Choi',ytitle='sps_fit'
  oplot,!x.crange,!x.crange,linestyle=2

  ploterror,all.agechoi,agediff,all.agechoierr,agedifferr,xtitle='Choi',ytitle='sps_fit-Choi',psym=1
  params = linfit(all.agechoi,agediff,measure_errors=agedifferr,sigma=sigmaparams)
  oplot,!x.crange,[-2,-2],linestyle=1,color=fsc_color('cyan')
  oplot,!x.crange,[2,2],linestyle=1,color=fsc_color('cyan')
  ploterror,all.MgFechoi,agediff,all.MgFechoierr,agedifferr,xtitle='Choi [Mg/Fe]',ytitle='sps_fit-Choi',psym=1

  plothist,agediff,bin=0.5,ytitle='Galaxies',xtitle='Choi -sps_fit'
  xyouts,0.5,0.97,'AGE(Gyr), All wavelengths',charsize=1.3,alignment=0.5,/normal
  device,/close

  ;POLYCONT
   ;COMPARE THE MEASUREMENTS OF Z
  fehdiff = polycont.feh-polycont.fehchoi
  fehdifferr = sqrt(polycont.feherr^2+polycont.fehchoierr^2)
  set_plot,'ps'
  !p.multi=[0,2,2]
  !p.font=0
  psname = 'choi_metal_comparison_polycont.eps'
  device,filename=psname,xsize=20,ysize=15,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  ploterror,polycont.fehchoi,polycont.feh,polycont.fehchoierr,polycont.feherr,psym=1,xtitle='Choi',ytitle='sps_fit'
  oplot,!x.crange,!x.crange,linestyle=2

  ploterror,polycont.fehchoi,fehdiff,polycont.fehchoierr,fehdifferr,xtitle='Choi',ytitle='sps_fit-Choi ',psym=1
  params = linfit(polycont.fehchoi,fehdiff,measure_errors=fehdifferr,sigma=sigmaparams)
  oplot,!x.crange,[-0.1,-0.1],linestyle=1,color=fsc_color('cyan')
  oplot,!x.crange,[0.1,0.1],linestyle=1,color=fsc_color('cyan')

  ploterror,polycont.MgFechoi,fehdiff,polycont.MgFechoierr,fehdifferr,xtitle='Choi [Mg/Fe]',ytitle='sps_fit-Choi ',psym=1

  plothist,fehdiff,bin=0.05,ytitle='Galaxies',xtitle='Choi-sps_fit'
  xyouts,0.5,0.97,'[Fe/H], Polynomial Normalization',charsize=1.3,alignment=0.5,/normal
  device,/close

;COMPARE THE MEASUREMENTS OF AGE
  agediff = polycont.age-polycont.agechoi
  agedifferr = sqrt(polycont.ageerr^2+polycont.agechoierr^2)

  !p.multi=[0,2,2]
  !p.font=0
  psname = 'choi_age_comparison_polycont.eps'
  device,filename=psname,xsize=20,ysize=15,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  ploterror,polycont.agechoi,polycont.age,polycont.agechoierr,polycont.ageerr,psym=1,xtitle='Choi',ytitle='sps_fit'
  oplot,!x.crange,!x.crange,linestyle=2

  ploterror,polycont.agechoi,agediff,polycont.agechoierr,agedifferr,xtitle='Choi',ytitle='sps_fit-Choi',psym=1
  params = linfit(polycont.agechoi,agediff,measure_errors=agedifferr,sigma=sigmaparams)
  oplot,!x.crange,[-2,-2],linestyle=1,color=fsc_color('cyan')
  oplot,!x.crange,[2,2],linestyle=1,color=fsc_color('cyan')
  ploterror,polycont.fehchoi,agediff,polycont.fehchoierr,agedifferr,xtitle='Choi [Fe/H]',ytitle='sps_fit-Choi',psym=1

  plothist,agediff,bin=0.5,ytitle='Galaxies',xtitle='Choi -sps_fit'
  xyouts,0.5,0.97,'AGE(Gyr), Polynomial Normalization',charsize=1.3,alignment=0.5,/normal
  device,/close

;POLYCONTV
   ;COMPARE THE MEASUREMENTS OF Z
  fehdiff = polycontv.feh-polycontv.fehchoi
  fehdifferr = sqrt(polycontv.feherr^2+polycontv.fehchoierr^2)
  set_plot,'ps'
  !p.multi=[0,2,2]
  !p.font=0
  psname = 'choi_metal_comparison_polycontv.eps'
  device,filename=psname,xsize=20,ysize=15,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  ploterror,polycontv.fehchoi,polycontv.feh,polycontv.fehchoierr,polycontv.feherr,psym=1,xtitle='Choi',ytitle='sps_fit'
  oplot,!x.crange,!x.crange,linestyle=2

  ploterror,polycontv.fehchoi,fehdiff,polycontv.fehchoierr,fehdifferr,xtitle='Choi',ytitle='sps_fit-Choi ',psym=1
  params = linfit(polycontv.fehchoi,fehdiff,measure_errors=fehdifferr,sigma=sigmaparams)
  ;oplot,!x.crange,params[0]+params[1]*!x.crange,color=fsc_color('darkgreen')
  oplot,!x.crange,[-0.1,-0.1],linestyle=1,color=fsc_color('cyan')
  oplot,!x.crange,[0.1,0.1],linestyle=1,color=fsc_color('cyan')

  ploterror,polycontv.MgFechoi,fehdiff,polycontv.MgFechoierr,fehdifferr,xtitle='Choi [Mg/Fe]',ytitle='sps_fit-Choi ',psym=1

  plothist,fehdiff,bin=0.05,ytitle='Galaxies',xtitle='Choi-sps_fit'
  xyouts,0.5,0.97,'[Fe/H], Polynomial Normalization (Vary n)',charsize=1.3,alignment=0.5,/normal
  device,/close

;COMPARE THE MEASUREMENTS OF AGE
  agediff = polycontv.age-polycontv.agechoi
  agedifferr = sqrt(polycontv.ageerr^2+polycontv.agechoierr^2)

  !p.multi=[0,2,2]
  !p.font=0
  psname = 'choi_age_comparison_polycontv.eps'
  device,filename=psname,xsize=20,ysize=15,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  ploterror,polycontv.agechoi,polycontv.age,polycontv.agechoierr,polycontv.ageerr,psym=1,xtitle='Choi',ytitle='sps_fit'
  oplot,!x.crange,!x.crange,linestyle=2

  ploterror,polycontv.agechoi,agediff,polycontv.agechoierr,agedifferr,xtitle='Choi',ytitle='sps_fit-Choi',psym=1
  params = linfit(polycontv.agechoi,agediff,measure_errors=agedifferr,sigma=sigmaparams)
  ;oplot,!x.crange,params[0]+params[1]*!x.crange,color=fsc_color('darkgreen')
  oplot,!x.crange,[-2,-2],linestyle=1,color=fsc_color('cyan')
  oplot,!x.crange,[2,2],linestyle=1,color=fsc_color('cyan')
  ploterror,polycontv.fehchoi,agediff,polycontv.fehchoierr,agedifferr,xtitle='Choi [Fe/H]',ytitle='sps_fit-Choi',psym=1

  plothist,agediff,bin=0.5,ytitle='Galaxies',xtitle='Choi -sps_fit'
  xyouts,0.5,0.97,'AGE(Gyr), Polynomial Normalization (Vary n)',charsize=1.3,alignment=0.5,/normal
  device,/close

  ;NO MGB
  ;COMPARE THE MEASUREMENTS OF Z
  fehdiff = noMgb.feh-noMgb.fehchoi
  fehdifferr = sqrt(noMgb.feherr^2+noMgb.fehchoierr^2)

  !p.multi=[0,2,2]
  !p.font=0
  psname = 'choi_metal_comparison_noMgb.eps'
  device,filename=psname,xsize=20,ysize=15,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  ploterror,noMgb.fehchoi,noMgb.feh,noMgb.fehchoierr,noMgb.feherr,psym=1,xtitle='Choi ',ytitle='sps_fit  (noMgb wl)'
  oplot,!x.crange,!x.crange,linestyle=2

  ploterror,noMgb.fehchoi,fehdiff,noMgb.fehchoierr,fehdifferr,xtitle='Choi ',ytitle='sps_fit (NoMgb wl) -Choi ',psym=1
  params = linfit(noMgb.fehchoi,fehdiff,measure_errors=fehdifferr,sigma=sigmaparams)
  ;oplot,!x.crange,params[0]+params[1]*!x.crange,color=fsc_color('darkgreen')
  oplot,!x.crange,[-0.1,-0.1],linestyle=1,color=fsc_color('cyan')
  oplot,!x.crange,[0.1,0.1],linestyle=1,color=fsc_color('cyan')
  ploterror,noMgb.MgFechoi,fehdiff,noMgb.MgFechoierr,fehdifferr,xtitle='Choi [Mg/Fe]',ytitle='sps_fit (noMgb wl) -Choi ',psym=1

  plothist,fehdiff,bin=0.05,ytitle='Galaxies',xtitle='Choi-sps_fit(noMgb)'
  device,/close

;COMPARE THE MEASUREMENTS OF AGE
  agediff = noMgb.age-noMgb.agechoi
  agedifferr = sqrt(noMgb.ageerr^2+noMgb.agechoierr^2)

  !p.multi=[0,2,2]
  !p.font=0
  psname = 'choi_age_comparison_noMgb.eps'
  device,filename=psname,xsize=20,ysize=15,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  ploterror,noMgb.agechoi,noMgb.age,noMgb.agechoierr,noMgb.ageerr,psym=1,xtitle='Choi Age',ytitle='sps_fit Age (noMgb wl)'
  oplot,!x.crange,!x.crange,linestyle=2

  ploterror,noMgb.agechoi,agediff,noMgb.agechoierr,agedifferr,xtitle='Choi Age',ytitle='sps_fit (NoMgb wl) Age-Choi Age',psym=1
  params = linfit(noMgb.agechoi,agediff,measure_errors=agedifferr,sigma=sigmaparams)
  ;oplot,!x.crange,params[0]+params[1]*!x.crange,color=fsc_color('darkgreen')
  oplot,!x.crange,[-2,-2],linestyle=1,color=fsc_color('cyan')
  oplot,!x.crange,[2,2],linestyle=1,color=fsc_color('cyan')
  ploterror,noMgb.MgFechoi,agediff,noMgb.MgFechoierr,agedifferr,xtitle='Choi [Mg/Fe]',ytitle='sps_fit (noMgb wl) Age-Choi Age',psym=1

  plothist,agediff,bin=0.5,ytitle='Galaxies',xtitle='Choi Age -sps_fit(noMgb) Age'
  device,/close
  stop
end
