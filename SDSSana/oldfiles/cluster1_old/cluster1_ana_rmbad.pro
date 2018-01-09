pro cluster1_ana_rmbad
;For cluster 1 the redshift is around 0.102
;read science
  dir = '/scr2/nichal/workspace2/sps_fit/data/'
  scifits = dir+'cluster1/sps_fit.fits.gz'
  clustername = 'cluster1'
  scienceall = mrdfits(scifits,1,/silent)
;make histograms of all objects in the field
  zfit = scienceall.zfit
  zspec = scienceall.zspec
  cgood = where(zfit ne 0)
  psname = 'redshift_histogram.eps'
  set_plot,'ps'
  device, filename = psname,xsize = 10,ysize = 16, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  !p.multi=[0,1,3]
  !p.font=0
  plothist,zspec,xhist_zspec,yhist_zspec,bin=0.05,title='All galaxies within R200 with spectra'
  plothist,zspec(cgood),bin=0.01,title='Galaxies fitted with sps_fit (S/N)>6'
  plothist,zfit(cgood),bin=0.01,/overplot,color=fsc_color('red'),linestyle=2
  al_Legend,['z_sdss','z_sps_fit'],psym=[0,0],linestyle=[0,2],color=[fsc_color('black'),fsc_color('red')],box=0,/right,/top
  plothist,zspec(cgood),bin=0.005,xrange=[0.05,0.15],title='(Zoomed) Galaxies fitted with sps_fit (S/N)>6'
  plothist,zfit(cgood),bin=0.005,/overplot,color=fsc_color('red'),linestyle=2
  device,/close

;compare zfit and zspec
  psname = 'redshift_comparison.eps'
  set_plot,'ps'
  device, filename = psname,xsize = 15,ysize = 8, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  !p.multi=[0,1,1]
  good = where(zfit ne 0)
  plot,zspec(good),zfit(good)-zspec(good),xtitle='SDSS z',ytitle='sps_fit z-SDSS z',psym=1,title=clustername
  device,/close

;metallicity and age. compare with local galaxies.
  ;define member as those with 0.08<z<0.12
  member = where(zfit gt 0.08 and zfit lt 0.12,cmember)
  color='navy blue'
  if cmember gt 0 then begin
     feh = scienceall(member).feh
     feherr = scienceall(member).feherr
     age = scienceall(member).age
     ageerr = scienceall(member).ageerr
     lmass  = scienceall(member).logmstar
  endif
  goodclump1 = where(feh gt .538*lmass-6.311 and lmass lt 10.8)
  goodclump2 = where(feh gt -0.6 and lmass gt 10.8)
  goodclump = [goodclump1,goodclump2]
  goodmember = member(goodclump)
  feh = scienceall(goodmember).feh
  feherr = scienceall(goodmember).feherr
  age = scienceall(goodmember).age
  ageerr = scienceall(goodmember).ageerr
  lmass  = scienceall(goodmember).logmstar

  ;find massbin average
  psname = 'stellarMZR_rmbad.eps'
  set_plot,'ps'
  device, filename = psname,xsize = 15,ysize = 8, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  !p.multi=[0,1,1]
  plot,lmass,feh,psym=1,xtitle='log M',ytitle='[Fe/H]',title=clustername,xstyle=5,ystyle=5,/nodata

  massarr = [9.9,10.3,10.5,10.7,10.9,11.1,11.3,11.5]
  nmass = n_elements(massarr)-1
  marr = fltarr(nmass)+1./0
  marr_lo = fltarr(nmass)+1./0
  marr_hi = fltarr(nmass)+1./0
  feharr  = fltarr(nmass)+1./0
  fehstdarr = fltarr(nmass)+1./0
  feherrarr = fltarr(nmass)+1./0
  for mi=0,nmass-1 do begin
     marr(mi) = (massarr[mi]+massarr[mi+1])*0.5
     msel = where(lmass gt massarr[mi] and lmass lt massarr[mi+1],cmsel)
     print, cmsel
     if cmsel gt 1 then begin
        meanerr,feh(msel),feherr(msel),fehmean,sigmam,sigmad,sigmas
        feharr(mi) = fehmean
        feherrarr(mi) = sigmam
        fehstdarr(mi) = sigmas
        meanerr,lmass(msel),feherr(msel),mmean,a,b,c
        marr[mi] = mmean
        marr_lo[mi] = massarr[mi]
        marr_hi[mi] = massarr[mi+1]
     endif         
  endfor
  top = interpol(feharr+fehstdarr,marr,massarr)
  bottom = interpol(feharr-fehstdarr,marr,massarr)
  toswitch = where(bottom gt top,csw)
  if csw gt 0 then begin
     temp = top(toswitch)
     top(toswitch) = bottom(toswitch)
     bottom(toswitch) = temp
  endif
  x=[massarr,reverse(massarr)]
  y=[top,reverse(bottom)]
  polyfill,x,y,color=fsc_color('light yellow')
  cgerrplot,marr,feharr-feherrarr,feharr+feherrarr,color=fsc_color(color)
  cgerrplot,feharr,marr_lo,marr_hi,color=fsc_color(color),/horizontal

  plot,lmass,feh,psym=1,xtitle='log M',ytitle='[Fe/H]',title=clustername,xstyle=1,ystyle=1,/nodata,/noerase

  oploterror,lmass,feh,feherr,psym=1
  cgplot,lmass,feh,color=fsc_color(color),psym=16,symsize=0.5,/overplot

  cgplot,marr,feharr,color=fsc_color('green'),psym=14,symsize=1.5,/overplot
  device,/close

  ;AGE
  ;readage from FAST
  fastfile = '/scr2/nichal/workspace2/SDSSana/FAST/cluster1.fout'
  readcol,fastfile,id,z,ltau,metal,lage,Av,lmass,lsfr,lssfr,la2t,chi2,format='I,F,F,F,F,F,F,F,F,F,F'
  restore,'/scr2/nichal/workspace2/SDSSana/FAST/goodphot.sav'
  agefast = scienceall.age*1./0.
  agefast(goodphot) = 10.^(lage-9.) ;Gyr
  agefast = agefast(member)
  psname = 'Age_comparison_rmbad.eps'
  set_plot,'ps'
  device, filename = psname,xsize = 10,ysize = 14, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  !p.multi=[0,1,2]
  !p.font=0
  plot,agefast,age,psym=4,xtitle='Age from FAST (GYR)',ytitle='Age from sps_fit'
  plot,agefast,agefast-age,psym=4,xtitle='AGE from FAST',ytitle='Age(FAST-SPSFIT)'
  device,/close

end
