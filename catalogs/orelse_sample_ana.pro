pro orelse_sample_ana

  cat = mrdfits('Cl1604.fits.gz',1,hdr)
  dfile = file_search('/scr2/nichal/keck/deimos/Cl1604/deimos/spec1d*.fits.gz',count=cdfile)
  lfile =file_search('/scr2/nichal/keck/deimos/Cl1604/lris/spec1d*.fits.gz',count=clfile)
  mfile = file_search('/scr2/nichal/workspace2/sps_fit/data/Cl1604*/sps_fit.fits.gz',count=cmfile)
  
  
;1) get the matched catalog
  catmask = strtrim(cat.mask,2)
  catslit = strtrim(cat.slit,2)
  catobjname = strtrim(cat.phot_id,2)

  dfilebase = file_basename(dfile)
  dmatchpos = intarr(cdfile)
  for i=0,cdfile-1 do begin
     extensions = strsplit(dfilebase[i],'.',/extract)
     slit = extensions[2]
     mask = extensions[1]
     dmatchpos[i] = where(catmask eq mask and catslit eq slit,cm)
  endfor
  dcat = cat(dmatchpos)
  
  lfilebase = file_basename(lfile)
  lmatchpos = intarr(clfile)
  for i=0,clfile-1 do begin
     extensions = strsplit(lfilebase[i],'.',/extract)
     if n_elements(extensions) eq 6 then slit = extensions[2]
     if n_elements(extensions) eq 7 then slit = extensions[2]+'.'+extensions[3]
     mask = extensions[1]
     if strmid(mask,0,3) ne 'old' then mask='LRIS.'+mask
     if mask eq 'oldLRISd' then mask = 'oldLRIS'
     if strmid(slit,0,1) eq '0' then slit =  strmid(slit,1)
     lmatchpos[i] = where(catmask eq mask and catslit eq slit,cm)
  endfor
  badmatch = where(lmatchpos eq -1,cbadmatch)
  if cbadmatch gt 0 then begin
   print, lfilebase(badmatch), 'Not in Catalog'
   remove, badmatch, lfilebase, lfile, lmatchpos
endif
  lcat = cat(lmatchpos)
  
  mstr = []
  for i=0,n_elements(mfile)-1 do begin 
     strnow  = mrdfits(mfile[i],1)
     mstr = [mstr,strnow]
  endfor
  match,strtrim(mstr.objname,2),catobjname,w1,w2
  mstrall = mstr
  mstr = mstr(w1)
  mcat = cat(w2)
  wlris = where(mstr.lrisflag eq 1)
  wdeimos = where(mstr.deimosflag eq 1)
;stop
;2)plot redshift distribution

  set_plot,'ps'
  psname = 'orelse_stat.eps'
  device, filename = psname,xsize = 35,ysize = 18, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  !p.font=0
  plothist,dcat.z,bin=0.02,xrange=[0.82,0.98],/fill,fcolor=fsc_color('rose'),position=[0.1,0.58,0.35,0.98],xstyle=5,ystyle=4
  plothist,lcat.z,bin=0.02,/overplot,/fill,fcolor=fsc_color('lightcyan')
  plothist,mcat(wdeimos).z,bin=0.02,/overplot,fcolor=fsc_color('red'),/fline,forient=45,/fill,color=fsc_color('red'),linestyle=1
  plothist,mcat(wlris).z,bin=0.02,/overplot,fcolor=fsc_color('blue'),/fline,forient=135,/fill,color=fsc_color('blue'),linestyle=1
  plothist,dcat.z,bin=0.02,xrange=[0.82,.98],xtitle='redshift',position=[0.1,0.58,0.35,0.98],xstyle=1,/noerase

  oplot,[0.886,0.886],!y.crange,color=fsc_color('darkgreen'),linestyle=2
  al_legend,['DEIMOS','LRIS','DEIMOS+MOSFIRE','DEIMOS+LRIS'],psym=15,color=fsc_Color(['rose','lightcyan','red','blue'])
  
;3) plot mass distribution
  binmass = 0.5
  massrange = [8.5,12.]
  gooddmass = where(dcat.logmstar_sed_pfw gt 0.)
  plothist,dcat(gooddmass).logmstar_sed_pfw,bin=binmass,xrange=massrange,/fill,fcolor=fsc_color('rose'),position=[0.4,0.58,0.65,0.98],xstyle=5,ystyle=4,/noerase
  goodlmass = where(lcat.logmstar_sed_pfw gt 0.)
  plothist,lcat(goodlmass).logmstar_sed_pfw,bin=binmass,/overplot,/fill,fcolor=fsc_color('lightcyan')
  plothist,mcat(wdeimos).logmstar_sed_pfw,bin=binmass,/overplot,fcolor=fsc_color('red'),/fline,forient=45,/fill,color=fsc_color('red'),linestyle=1
  plothist,mcat(wlris).logmstar_sed_pfw,bin=binmass,/overplot,fcolor=fsc_color('blue'),/fline,forient=135,/fill,color=fsc_color('blue'),linestyle=1
  plothist,dcat(gooddmass).logmstar_sed_pfw,bin=binmass,xrange=massrange,xtitle='log[M*]',position=[0.4,0.58,0.65,0.98],xstyle=1,/noerase


;3) plot zmag distribution
  binz= 1
  zrange = [16,25.]
  gooddmass = where(dcat.zmag gt 0.)
  plothist,dcat(gooddmass).zmag,bin=binz,xrange=zrange,/fill,fcolor=fsc_color('rose'),position=[0.7,0.58,0.95,0.98],xstyle=5,ystyle=4,/noerase
  goodlmass = where(lcat.zmag gt 0.)
  plothist,lcat(goodlmass).zmag,bin=binz,/overplot,/fill,fcolor=fsc_color('lightcyan')
  plothist,mcat(wdeimos).zmag,bin=binz,/overplot,fcolor=fsc_color('red'),/fline,forient=45,/fill,color=fsc_color('red'),linestyle=1
  plothist,mcat(wlris).zmag,bin=binz,/overplot,fcolor=fsc_color('blue'),/fline,forient=135,/fill,color=fsc_color('blue'),linestyle=1
  plothist,dcat(gooddmass).zmag,bin=binz,xrange=zrange,xtitle='z mag',position=[0.7,0.58,0.95,0.98],xstyle=1,/noerase

;3) Observed with MOSFIRE
;3.1)
  plothist,mcat.z,bin=0.02,/noerase,position=[0.1,0.08,0.35,0.48],xstyle=1,xtitle='redshift',/ylog,yrange=[0.1,100]
  gooddeimos = where(mstr.gooddeimos eq 1)
  goodmosfire = where(mstr.goodmosfire eq 1 and mcat.zmag le 22.)
  goodall = where(mstr.gooddeimos eq 1 and mstr.goodmosfire eq 1)
  
  plothist,mcat(goodmosfire).z,xhist,yhist,bin=0.02,/noplot,/overplot
  vsym,3,/fill
  oplot,xhist,yhist,psym=8,color=fsc_color('orange')

  plothist,mcat(goodall).z,xhist,yhist,bin=0.02,/noplot,/overplot
  vsym,24,/fill
  oplot,xhist,yhist,psym=8,color=fsc_color('green')

  plothist,mcat(gooddeimos).z,xhist,yhist,bin=0.02,/noplot,/overplot
  vsym,5,/star,/fill
  oplot,xhist,yhist,psym=8,color=fsc_color('red')
  al_legend,['All MOSFIRE','Good DEIMOS/LRIS','Good MOSFIRE'],psym=[15,46,17],color=fsc_Color(['black','red','orange'])
  print, 'Good all spec', n_elements(goodall)
  
  ;3.2)
  plothist,mcat.logmstar_sed_pfw,bin=0.5,/noerase,position=[0.4,0.08,0.65,0.48],xstyle=1,xtitle='log [M*]',/ylog,yrange=[0.1,100]
  
  plothist,mcat(goodmosfire).logmstar_sed_pfw,xhist,yhist,bin=0.5,/noplot,/overplot
  vsym,3,/fill
  oplot,xhist,yhist,psym=8,color=fsc_color('orange')

  plothist,mcat(goodall).logmstar_sed_pfw,xhist,yhist,bin=0.5,/noplot,/overplot
  vsym,24,/fill
  oplot,xhist,yhist,psym=8,color=fsc_color('green')

  plothist,mcat(gooddeimos).logmstar_sed_pfw,xhist,yhist,bin=0.5,/noplot,/overplot
  vsym,5,/star,/fill
  oplot,xhist,yhist,psym=8,color=fsc_color('red')
  
;3.3)
  plothist,mcat.zmag,bin=1,/noerase,position=[0.7,0.08,0.95,0.48],xstyle=1,xtitle='log [M*]',/ylog,yrange=[0.1,100]
  
  plothist,mcat(goodmosfire).zmag,xhist,yhist,bin=1,/noplot,/overplot
  vsym,3,/fill
  oplot,xhist,yhist,psym=8,color=fsc_color('orange')

  plothist,mcat(goodall).zmag,xhist,yhist,bin=1,/noplot,/overplot
  vsym,24,/fill
  oplot,xhist,yhist,psym=8,color=fsc_color('green')

  plothist,mcat(gooddeimos).zmag,xhist,yhist,bin=1,/noplot,/overplot
  vsym,5,/star,/fill
  oplot,xhist,yhist,psym=8,color=fsc_color('red')

  device,/close
  
;PLOT Z MAG VS MASS
  psname = 'mass_zmag.eps'
  !p.multi=[0,1,1]
  device, filename = psname,xsize = 12,ysize = 10, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  plot,mcat.zmag,mcat.logmstar_sed_pfw,psym=1,xtitle='z mag', ytitle='log[M!D*!N/M!Dsun!N]',/nodata,yrange=[8.5,12]
  vsym,24,/fill
  oplot,mcat.zmag,mcat.logmstar_sed_pfw,color=fsc_color('darkgray'),psym=8
  vsym,3,/fill
  oplot,mcat(goodmosfire).zmag,mcat(goodmosfire).logmstar_sed_pfw,psym=8,color=fsc_color('orange')
  vsym,24,/fill
  oplot,mcat(goodall).zmag,mcat(goodall).logmstar_sed_pfw,psym=8,color=fsc_color('green')
  vsym,5,/star,/fill
  oplot,mcat(gooddeimos).zmag,mcat(gooddeimos).logmstar_sed_pfw,psym=8,color=fsc_color('red')
  vsym,5,/star,/fill
  al_legend,['All MOSFIRE&DEIMOS/LRIS', 'Good All','Good MOSFIRE','Good DEIMOS/LRIS'],psym=[16,16,17,8],color=fsc_color(['darkgray','green','orange','red']),/left,/bottom,box=0
  device,/close
  psname = 'mass_zmag2.eps'
  !p.multi=[0,1,1]
  device, filename = psname,xsize = 12,ysize = 10, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  plot,dcat(gooddmass).zmag,dcat(gooddmass).logmstar_sed_pfw,psym=1,xtitle='z mag', ytitle='log[M!D*!N/M!Isun!N]',/nodata,yrange=[8.5,12],xrange=[19,25]
  vsym,24,/fill
  oplot,dcat(gooddmass).zmag,dcat(gooddmass).logmstar_sed_pfw,psym=8,color=fsc_color('indianred')
  oplot,lcat(goodlmass).zmag,lcat(goodlmass).logmstar_sed_pfw,psym=8,color=fsc_color('cyan')
  al_legend,['All DEIMOS','All LRIS'],color=fsc_color(['indianred','cyan']),psym=16,/left,/bottom,box=0
  device,/close

;PLOT RA AND DEC
  ra0 = 241.1d
  dec0 = 43.30d
  rarange = [8, -8]
  decrange = [-8, 8]

  mosfirex = [1.5, 1.5, -1.5, -1.5, 1.5]
  mosfirey = [3.06, -3.06, -3.06, 3.06, 3.06]
  LRISx = [3.,3.,-3.,-3,3.]
  LRISy = [3.9,-3.9,-3.9,3.9,3.9]
  masks = 'Cl1604_'+['B','C','D1', 'D2']
  nmasks = n_elements(masks)
  raff = dblarr(nmasks)
  decff = dblarr(nmasks)
  paff = dblarr(nmasks)
  radec = ' '
 
  for i=0,nmasks-1 do begin
     openr, lun, masks[i]+'_StarList.txt', /get_lun
     readf, lun, radec, pa, format='(16X,A24,17X,D6)'
     close, lun
     free_lun, lun
     get_coords, coords, instring=radec
     raff[i] = coords[0]*15.
     decff[i] = coords[1]
     paff[i] = pa
  endfor
  nmasks=2
  decl = [43.252072,43.362667]
  ral = [241.04485,241.15517]
  pal = [-30.,67.5]
  decf = [43.235072,43.362667]
  raf = [241.09485,241.15517]
  paf = [5.,67.5]
  bright = where(mcat.zmag lt 22.)
  !p.multi = [0,1,1]
  device, filename='MOSFIRE_radec.eps', xsize=7, ysize=7, /inches, /color, /encapsulated
  vsym,24,/fill  
  plot, (mcat.ra-ra0)*cos(dec0*!DTOR)*60., (mcat.dec-dec0)*60., psym=8, xtitle='!Ma !M-!7 !Ma!7!D0!N', ytitle='!Md !M-!7 !Md!7!D0!N', xrange=rarange, yrange=decrange, xstyle=1, ystyle=1, /isotropic,/nodata
  oplot,(mcat.ra-ra0)*cos(dec0*!DTOR)*60., (mcat.dec-dec0)*60., psym=8,color=fsc_color('darkgray')
  oplot,(mcat(bright).ra-ra0)*cos(dec0*!DTOR)*60., (mcat(bright).dec-dec0)*60., psym=8,color=fsc_color('darkgray')

  vsym,3,/fill
  oplot,(mcat(goodmosfire).ra-ra0)*cos(dec0*!DTOR)*60.,(mcat(goodmosfire).dec-dec0)*60.,psym=8,color=fsc_color('orange')
  vsym,24,/fill
  oplot,(mcat(goodall).ra-ra0)*cos(dec0*!DTOR)*60.,(mcat(goodall).dec-dec0)*60.,psym=8,color=fsc_color('green')
  vsym,5,/star,/fill
  oplot,(mcat(gooddeimos).ra-ra0)*cos(dec0*!DTOR)*60.,(mcat(gooddeimos).dec-dec0)*60.,psym=8,color=fsc_color('red')
  vsym,5,/star,/fill

  for i=0,nmasks-1 do oplot, mosfirex*cos(paf[i]*!DTOR)+mosfirey*sin(paf[i]*!DTOR)+(raf[i]-ra0)*cos(dec0*!DTOR)*60., -1.0*mosfirex*sin(paf[i]*!DTOR)+mosfirey*cos(paf[i]*!DTOR)+(decf[i]-dec0)*60., color=fsc_color('blue'), thick=8
  for i=0,nmasks-1 do oplot, lrisx*cos(pal[i]*!DTOR)+lrisy*sin(pal[i]*!DTOR)+(ral[i]-ra0)*cos(dec0*!DTOR)*60., -1.0*lrisx*sin(pal[i]*!DTOR)+lrisy*cos(pal[i]*!DTOR)+(decl[i]-dec0)*60., color=fsc_color('orange'), thick=8
  xyouts,1,7,'A',color=fsc_color('orange'),charsize=1.5
  xyouts,2,5.5,'C',color=fsc_color('blue'),charsize=1.5
  xyouts,-5,.5,'B',color=fsc_color('orange'),charsize=1.5
  xyouts,2.5,-3,'D',color=fsc_color('blue'),charsize=1.5

  al_legend,['All MOSFIRE&DEIMOS/LRIS', 'Good All','Good MOSFIRE','Good DEIMOS/LRIS'],psym=[16,16,17,8],color=fsc_color(['darkgray','green','orange','red']),/left,/bottom,box=0
  device,/close


  device, filename='MOSFIRE_radec_2018a.eps', xsize=7, ysize=7, /inches, /color, /encapsulated
  vsym,24,/fill
  plot, (mcat.ra-ra0)*cos(dec0*!DTOR)*60., (mcat.dec-dec0)*60., psym=8, xtitle='!Ma !M-!7 !Ma!7!D0!N', ytitle='!Md !M-!7 !Md!7!D0!N', xrange=rarange, yrange=decrange, xstyle=1, ystyle=1, /isotropic,/nodata
  oplot,(mcat.ra-ra0)*cos(dec0*!DTOR)*60., (mcat.dec-dec0)*60., psym=8,color=fsc_color('darkgray')
  oplot,(mcat(bright).ra-ra0)*cos(dec0*!DTOR)*60., (mcat(bright).dec-dec0)*60., psym=8,color=fsc_color('darkgray')

  vsym,3,/fill
  oplot,(mcat(goodmosfire).ra-ra0)*cos(dec0*!DTOR)*60.,(mcat(goodmosfire).dec-dec0)*60.,psym=8,color=fsc_color('orange')
  vsym,24,/fill
  oplot,(mcat(goodall).ra-ra0)*cos(dec0*!DTOR)*60.,(mcat(goodall).dec-dec0)*60.,psym=8,color=fsc_color('green')
  vsym,5,/star,/fill
  oplot,(mcat(gooddeimos).ra-ra0)*cos(dec0*!DTOR)*60.,(mcat(gooddeimos).dec-dec0)*60.,psym=8,color=fsc_color('red')
  vsym,5,/star,/fill

  for i=1,nmasks-1 do oplot, mosfirex*cos(paf[i]*!DTOR)+mosfirey*sin(paf[i]*!DTOR)+(raf[i]-ra0)*cos(dec0*!DTOR)*60., -1.0*mosfirex*sin(paf[i]*!DTOR)+mosfirey*cos(paf[i]*!DTOR)+(decf[i]-dec0)*60., color=fsc_color('blue'), thick=8
  for i=0,nmasks-1 do oplot, lrisx*cos(pal[i]*!DTOR)+lrisy*sin(pal[i]*!DTOR)+(ral[i]-ra0)*cos(dec0*!DTOR)*60., -1.0*lrisx*sin(pal[i]*!DTOR)+lrisy*cos(pal[i]*!DTOR)+(decl[i]-dec0)*60., color=fsc_color('orange'), thick=8
  xyouts,1,7,'A',color=fsc_color('orange'),charsize=1.5
  xyouts,2,5.5,'C',color=fsc_color('blue'),charsize=1.5
  xyouts,-5,.5,'B',color=fsc_color('orange'),charsize=1.5
;  xyouts,2.5,-3,'D',color=fsc_color('blue'),charsize=1.5

  al_legend,['All MOSFIRE&DEIMOS/LRIS', 'Good All','Good MOSFIRE','Good DEIMOS/LRIS'],psym=[16,16,17,8],color=fsc_color(['darkgray','green','orange','red']),/left,/bottom,box=0

  device,/close


  stop
end
