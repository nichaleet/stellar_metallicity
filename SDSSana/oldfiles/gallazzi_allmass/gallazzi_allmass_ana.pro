pro gallazzi_allmass_ana
  ;READ DATA
  ;read science
  sci = mrdfits('/scr2/nichal/workspace2/sps_fit/data/gallazzi_allmass/sps_fit.fits.gz',1,/silent)
	sci.age = alog10(sci.age)+9.
	sci.ageerr = sci.ageerr/sci.age/alog(10.)
  ;read catalog
	cat = mrdfits('/scr2/nichal/workspace2/SDSSdata/gallazzi_allmass/z_mstar_morphology.fits',1)
	;get the match
	nobjs = n_elements(Sci)
  matcharr = lonarr(nobjs)
  for i=0,nobjs-1 do begin
     match = where(sci[i].plate eq cat.plate and sci[i].mjd eq cat.mjd and sci[i].fiber eq cat.fiber,cmatch)
     if cmatch ne 1 then stop
     matcharr[i] = match
  endfor
	cat = cat(matcharr)
	;check the match
	if total(cat.plate-sci.plate) ne 0 or total(cat.fiber-cat.fiber) ne 0 then stop,'check the match'

  ;fix the age see Gallazzi06
 ; added_age = (galage(0,1000.)-galage(sci.zfit,1000.))/1.e9
 ; cat.agelow = cat.agelow-added_age
 ; cat.age = cat.age-added_age
 ; cat.agehigh = cat.agehigh-added_age
  ;stop

  save,cat,filename='match_gallazzi_allmass.sav'
  ;stop

  ;REDSHIFT,MASS DISTRIBUTION
  set_plot,'ps'
  !p.multi = [0,2,1]
  ;!p.font = 0
  psname = 'redshift_mass_hist_allmass.eps'
  device, filename = psname,xsize = 20,ysize = 8, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  plothist,sci.zfit,bin=0.01,xtitle='Redshift',xtickformat='(F3.1)'
  plothist,cat.logm50,bin=0.1,xtitle='Log(Mass)',xrange=[9,13]

  device,/close

  ;COMPARE THE MEASUREMENTS OF Z and AGE
  fehdiff = sci.feh-cat.z50
  ;fix sigma
  badfeherr = where(sci.feherr eq 0., cbadfeherr)
  if cbadfeherr gt 0 then sci(badfeherr).feherr = median(sci.feherr)

  fehdiff_sigma = (sci.feh-cat.z50)/sqrt(sci.feherr^2+(0.5*(cat.z84-cat.z16))^2);number of sigmas differences
  fehdifflow = fehdiff-sqrt(sci.feherr^2+(cat.z50-cat.z16)^2)
  fehdiffhigh = fehdiff+sqrt(sci.feherr^2+(cat.z50-cat.z84)^2)

  agediff = sci.age-cat.logt50
  agediff_sigma = (sci.age-cat.logt50)/sqrt(sci.ageerr^2+(0.5*(cat.logt84-cat.logt16))^2)
  agedifflow = agediff-sqrt(sci.ageerr^2+(cat.logt50-cat.logt16)^2)
  agediffhigh = agediff+sqrt(sci.ageerr^2+(cat.logt50-cat.logt84)^2)
 
  noem = where(sci.haveemlines eq 0, cnoem) ;no emission lines
  em = where(sci.haveemlines eq 1, cem) ; with emission lines

  !p.multi=[0,1,1]
  !p.charsize=1
  psname = 'age_metal_comparison_sigma_allmass.eps'
  device,filename=psname,xsize=35,ysize=10,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  ;1) FEH diff (sigma)
  Delta = '!4'+string("104B)+'!x'
  sigma = '!4'+string("162B)+'!x'
  plothist,fehdiff_sigma,xtitle=Delta+' [Fe/H]('+sigma+')',position=[0.1,0.1,0.33,0.95],bin=1,/halfbin
  plothist,fehdiff_sigma(em),/overplot,color=fsc_color('darkgreen'),/fill,fcolor=fsc_color('cyan'),bin=1,/halfbin
  plothist,fehdiff_sigma(noem),/overplot,color=fsc_color('red'),bin=1,/halfbin

  ;2) AGE diff (sigma)
  plothist,agediff_sigma,xtitle=Delta+'log(Age/yr) ('+sigma+')',position=[0.43,0.1,0.66,0.95],/noerase,bin=1,/halfbin
  plothist,agediff_sigma(em),/overplot,color=fsc_color('darkgreen'),/fill,fcolor=fsc_color('cyan'),bin=1,/halfbin
  plothist,agediff_sigma(noem),/overplot,color=fsc_color('red'),bin=1,/halfbin

  ;3) AGE vs FEH
  ;;linear fit
  linparam = linfit(agediff,fehdiff)
  print,'slope of the correlation:',linparam(1)
  plot,fehdiff,agediff,psym=1,xtitle=delta+'[Fe/H](dex)',ytitle=delta+'log(Age)',/nodata,/noerase,position=[0.76,0.1,0.97,0.95]
  cgerrplot,agediff,fehdifflow,fehdiffhigh,/horizontal,color='darkgray'
  cgerrplot,fehdiff,agedifflow,agediffhigh,color='darkgray'
  vsym,4,/fill,rot=45
  oplot,fehdiff(noem),agediff(noem),psym=8,color=fsc_color('pink')
  oplot,fehdiff(em),agediff(em),psym=8,color=fsc_color('cyan')
  al_legend,['no Emission lines','with Emission lines'],psym=15,color=fsc_color(['pink','cyan'])

  device,/close    	
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
  !P.MULTI = 0
  !p.font=0
  !p.charsize=0
  !p.position = [0.1,0.15,0.95,0.95]
  Delta = '!9'+string("104B)+'!x'
  Delta = '!9'+string("104B)+'!x'

  psname = 'feh_nl_gallazzi_allmass.eps'
  device,filename=psname,xsize=12,ysize=8,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color

  plothist,fehdiff,bin=0.1,/halfbin,/fill,fcolor=fsc_color('blk3'),xrange=[-0.8,0.8],xstyle=5,yrange=[0,60],ystyle=5
  plothist,fehdiff(noem),fehdiff_xhist,fehdiff_yhist,xtitle=Delta+' Age ('+sigma+')',/overplot,color=fsc_color('org6'),bin=0.1,/halfbin,/fill,fcolor=fsc_Color('org4'),thick=2
  plothist,fehdiff(em),/overplot,color=fsc_color('blu7'),bin=0.1,/halfbin,/fline,forient=45,fcolor=fsc_color('blu5'),thick=2
	gauss_noem = gaussfit(fehdiff_xhist,fehdiff_yhist,coeff,nterms=3)
	xfit = findgen(101)/50.-1.
	gauss_noem = gaussian(xfit,coeff)
	oplot,xfit,gauss_noem,color=fsc_color('org6')
 	axis,xaxis=0,xtitle=Delta+' [Fe/H](dex)'
	axis,xaxis=1,xtickformat='(A1)'
	axis,yaxis=1,ytickformat='(A1)'
	axis,yaxis=0
  al_legend,['all','passive','star forming'],psym=10,color=fsc_color(['black','org6','blu7']),position=[-0.8,55],box=0,/poly_fill,polycolor=['blk3','org4','blu5'],line_orientation=[-200,-200,45],polyspace=0.1
	device,/close

	psname = 'age_nl_gallazzi_allmass.eps'
  device,filename=psname,xsize=12,ysize=8,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   plothist,agediff,xtitle=delta+'log(Age)',bin=0.1,/halfbin,/fill,fcolor=('blk3'),xstyle=5,ystyle=4,xrange=[-0.5,0.4]
;  plothist,agediff/cat.age,xtitle=Delta+' Age (Gyr)/Aige',bin=0.1,/halfbin,/fill,fcolor=('blk3'),xstyle=4,ystyle=4
   plothist,agediff(noem),/overplot,color=fsc_color('org6'),bin=0.1,/halfbin,/fill,fcolor=fsc_color('org4'),thick=2
;  plothist,agediff(noem)/cat.age(noem),/overplot,color=fsc_color('org6'),bin=0.1,/halfbin,/fill,fcolor=fsc_color('org4'),thick=2
   plothist,agediff(em),/overplot,color=fsc_color('blu7'),/fline,forient=45,fcolor=fsc_color('blu5'),bin=0.1,/halfbin,thick=2
  ;plothist,agediff(em)/cat.age(em),/overplot,color=fsc_color('blu7'),/fline,forient=45,fcolor=fsc_color('blu5'),bin=0.1,/halfbin,thick=2
  axis,xaxis=0,xtitle=Delta+'log(Age)',xrange=[-0.5,0.4],xstyle=1
  axis,xaxis=1,xtickformat='(A1)'
  axis,yaxis=1,ytickformat='(A1)'
  axis,yaxis=0
	
 plot,sci(noem).feh,agediff(noem),psym=1,xtitle='[Fe/H](dex)',ytitle=delta+'log(Age)',/nodata,/noerase,position=[0.6,0.6,0.93,0.93],charsize=0.8
  cgerrplot,agediff(noem),sci(noem).feh-sci(noem).feherr,sci(noem).feh+sci(noem).feherr,/horizontal,color='org4'
  cgerrplot,sci(noem).feh,agedifflow(noem),agediffhigh(noem),color='org4'
  vsym,4,/fill,rot=45
  oplot,sci(noem).feh,agediff(noem),psym=8,color=fsc_color('org4')
	linpar = linfit(sci(noem).feh,agediff(noem),measure_errors=(agediffhigh(noem)-agedifflow(noem))/2.)
  ;oplot,!x.crange,linpar(0)+!x.crange*linpar(1),color=fsc_color('org6')

  device,/close

;  psname = 'metal_comparison_histogram.eps'
;	cleanplot
;	!p.multi=[0,3,1]
;  device,filename=psname,xsize=20,ysize=8,$
;         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
;	plot,cat.age,sci.age,/nodata,xrange=[0,14],xtitle='Gallazzi Age',ytitle='NL age'
;	oplot,cat.age(em),sci(em).age,color=fsc_color('blu5'),psym=1
;	oplot,cat.age(noem),sci(noem).age,color=fsc_color('org4'),psym=1
;	oplot,!x.crange,!x.crange
;  plot,cat.feh,sci.feh,/nodata,xtitle='Gallazzi Feh',xrange=[-0.5,0.5]
;  oplot,cat.feh(em),sci(em).feh,color=fsc_color('blu5'),psym=1
;  oplot,cat.feh(noem),sci(noem).feh,color=fsc_color('org4'),psym=1
;  oplot,!x.crange,!x.crange
;  plot,sci.feh,alog10(sci.age)-alog10(cat.age),/nodata,xtitle='FeH',ytitle='log(AGE_NL)-log(AGE_gal)'
;  oplot,sci(em).feh,alog10(sci(em).age)-alog10(cat.age(em)),color=fsc_Color('blu5'),psym=1
;  oplot,sci(noem).feh,alog10(sci(noem).age)-alog10(cat.age(noem)),color=fsc_Color('org4'),psym=1	
;  device,/close
stop
end
