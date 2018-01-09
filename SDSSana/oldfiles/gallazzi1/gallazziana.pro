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
  ;fix the age see Gallazzi06
 ; added_age = (galage(0,1000.)-galage(sci.zfit,1000.))/1.e9
 ; cat.agelow = cat.agelow-added_age
 ; cat.age = cat.age-added_age
 ; cat.agehigh = cat.agehigh-added_age
  ;stop
  save,cat,filename='match_gallazzi1.sav'
  ;stop

  ;REDSHIFT,MASS DISTRIBUTION
  set_plot,'ps'
  !p.multi = [0,2,1]
  ;!p.font = 0
  psname = 'redshift_mass_hist.eps'
  device, filename = psname,xsize = 20,ysize = 8, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  plothist,sci.zfit,bin=0.01,xtitle='Redshift',xtickformat='(F3.1)'
  plothist,cat.logmass,bin=0.1,xtitle='Log(Mass)',xrange=[9,13]

  device,/close

  ;COMPARE THE MEASUREMENTS OF Z and AGE
  fehdiff = sci.feh-cat.feh
  ;fix sigma
  badfeherr = where(sci.feherr eq 0., cbadfeherr)
  if cbadfeherr gt 0 then sci(badfeherr).feherr = median(sci.feherr)
  fehdiff_sigma = (sci.feh-cat.feh)/sqrt(sci.feherr^2+(0.5*(cat.fehhigh-cat.fehlow))^2)
  fehdifflow = fehdiff-sqrt(sci.feherr^2+(cat.feh-cat.fehlow)^2)
  fehdiffhigh = fehdiff+sqrt(sci.feherr^2+(cat.feh-cat.fehhigh)^2)

  agediff = sci.age-cat.age
  agediff_sigma = (sci.age-cat.age)/sqrt(sci.ageerr^2+(0.5*(cat.agehigh-cat.agelow))^2)
  agedifflow = agediff-sqrt(sci.ageerr^2+(cat.age-cat.agelow)^2)
  agediffhigh = agediff+sqrt(sci.ageerr^2+(cat.age-cat.agehigh)^2)
 
  noem = where(sci.haveemlines eq 0, cnoem) ;no emission lines
  em = where(sci.haveemlines eq 1, cem) ; with emission lines

  !p.multi=[0,1,1]
  !p.charsize=1
  psname = 'age_metal_comparison_sigma.eps'
  device,filename=psname,xsize=35,ysize=10,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  ;1) FEH diff (sigma)
  Delta = '!4'+string("104B)+'!x'
  sigma = '!4'+string("162B)+'!x'
  plothist,fehdiff_sigma,xtitle=Delta+' [Fe/H]('+sigma+')',position=[0.1,0.1,0.33,0.95],bin=1,/halfbin
  plothist,fehdiff_sigma(em),xtitle=Delta+' Age ('+sigma+')',/overplot,color=fsc_color('darkgreen'),/fill,fcolor=fsc_color('cyan'),bin=1,/halfbin
  plothist,fehdiff_sigma(noem),xtitle=Delta+' Age ('+sigma+')',/overplot,color=fsc_color('red'),bin=1,/halfbin

  ;2) AGE diff (sigma)
  plothist,agediff_sigma,xtitle=Delta+' Age ('+sigma+')',position=[0.43,0.1,0.66,0.95],/noerase,bin=1,/halfbin
  plothist,agediff_sigma(em),xtitle=Delta+' Age ('+sigma+')',/overplot,color=fsc_color('darkgreen'),/fill,fcolor=fsc_color('cyan'),bin=1,/halfbin
  plothist,agediff_sigma(noem),xtitle=Delta+' Age ('+sigma+')',/overplot,color=fsc_color('red'),bin=1,/halfbin

  ;3) AGE vs FEH
  logagediff = alog10(sci.age)-alog10(cat.age)
  logagediffhigh = logagediff+sqrt((sci.ageerr/sci.age/alog(10.))^2+((cat.agehigh-cat.age)/cat.age/alog(10.))^2)
  logagedifflow = logagediff-sqrt((sci.ageerr/sci.age/alog(10.))^2+((cat.age-cat.agelow)/cat.age/alog(10.))^2)
  ;;linear fit
  linparam = linfit(logagediff,fehdiff)
  print,'slope of the correlation:',linparam(1)
  plot,fehdiff,logagediff,psym=1,xtitle=delta+'[Fe/H](dex)',ytitle=delta+'log(Age)',/nodata,/noerase,position=[0.76,0.1,0.97,0.95]
  cgerrplot,logagediff,fehdifflow,fehdiffhigh,/horizontal,color='darkgray'
  cgerrplot,fehdiff,logagedifflow,logagediffhigh,color='darkgray'
  vsym,4,/fill,rot=45
  oplot,fehdiff(noem),logagediff(noem),psym=8,color=fsc_color('pink')
  oplot,fehdiff(em),logagediff(em),psym=8,color=fsc_color('cyan')
  al_legend,['no Emission lines','with Emission lines'],psym=15,color=fsc_color(['pink','cyan'])

  device,/close    	
	
  ;3) FEH diff (dex)
  !P.MULTI = 0
  !p.font=0
  !p.charsize=0
  !p.position = [0.1,0.15,0.95,0.95]
  Delta = '!9'+string("104B)+'!x'
  Delta = '!9'+string("104B)+'!x'

  psname = 'feh_nl_gallazzi.eps'
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

  ;4) AGE diff (Gyr)
	psname = 'age_nl_gallazzi.eps'
  device,filename=psname,xsize=12,ysize=8,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   plothist,alog10(sci.age)-alog10(cat.age),xtitle=delta+'log(Age)',bin=0.1,/halfbin,/fill,fcolor=('blk3'),xstyle=5,ystyle=4,xrange=[-0.5,0.4]
;  plothist,agediff/cat.age,xtitle=Delta+' Age (Gyr)/Aige',bin=0.1,/halfbin,/fill,fcolor=('blk3'),xstyle=4,ystyle=4
   plothist,alog10(sci(noem).age)-alog10(cat.age(noem)),/overplot,color=fsc_color('org6'),bin=0.1,/halfbin,/fill,fcolor=fsc_color('org4'),thick=2
;  plothist,agediff(noem)/cat.age(noem),/overplot,color=fsc_color('org6'),bin=0.1,/halfbin,/fill,fcolor=fsc_color('org4'),thick=2
   plothist,alog10(sci(em).age)-alog10(cat.age(em)),/overplot,color=fsc_color('blu7'),/fline,forient=45,fcolor=fsc_color('blu5'),bin=0.1,/halfbin,thick=2
  ;plothist,agediff(em)/cat.age(em),/overplot,color=fsc_color('blu7'),/fline,forient=45,fcolor=fsc_color('blu5'),bin=0.1,/halfbin,thick=2
  axis,xaxis=0,xtitle=Delta+'log(Age)',xrange=[-0.5,0.4],xstyle=1
  axis,xaxis=1,xtickformat='(A1)'
  axis,yaxis=1,ytickformat='(A1)'
  axis,yaxis=0
	
 plot,sci(noem).feh,logagediff(noem),psym=1,xtitle='[Fe/H](dex)',ytitle=delta+'log(Age)',/nodata,/noerase,position=[0.6,0.6,0.93,0.93],charsize=0.8
  cgerrplot,logagediff(noem),sci(noem).feh-sci(noem).feherr,sci(noem).feh+sci(noem).feherr,/horizontal,color='org4'
  cgerrplot,sci(noem).feh,logagedifflow(noem),logagediffhigh(noem),color='org4'
  vsym,4,/fill,rot=45
  oplot,sci(noem).feh,logagediff(noem),psym=8,color=fsc_color('org4')
	linpar = linfit(sci(noem).feh,logagediff(noem),measure_errors=(logagediffhigh(noem)-logagedifflow(noem))/2.)
  ;oplot,!x.crange,linpar(0)+!x.crange*linpar(1),color=fsc_color('org6')

  device,/close

  psname = 'metal_comparison_histogram.eps'
	cleanplot
	!p.multi=[0,3,1]
  device,filename=psname,xsize=20,ysize=8,$
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
	plot,cat.age,sci.age,/nodata,xrange=[0,14],xtitle='Gallazzi Age',ytitle='NL age'
	oplot,cat.age(em),sci(em).age,color=fsc_color('blu5'),psym=1
	oplot,cat.age(noem),sci(noem).age,color=fsc_color('org4'),psym=1
	oplot,!x.crange,!x.crange
  plot,cat.feh,sci.feh,/nodata,xtitle='Gallazzi Feh',xrange=[-0.5,0.5]
  oplot,cat.feh(em),sci(em).feh,color=fsc_color('blu5'),psym=1
  oplot,cat.feh(noem),sci(noem).feh,color=fsc_color('org4'),psym=1
  oplot,!x.crange,!x.crange
  plot,sci.feh,alog10(sci.age)-alog10(cat.age),/nodata,xtitle='FeH',ytitle='log(AGE_NL)-log(AGE_gal)'
  oplot,sci(em).feh,alog10(sci(em).age)-alog10(cat.age(em)),color=fsc_Color('blu5'),psym=1
  oplot,sci(noem).feh,alog10(sci(noem).age)-alog10(cat.age(noem)),color=fsc_Color('org4'),psym=1	
  device,/close
stop
end
