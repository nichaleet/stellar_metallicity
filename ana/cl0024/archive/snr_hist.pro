pro snr_hist,maskarr,zmin,zmax
;e.g. snr_hist,['rse1','rse2','rse3','rse4','rse5','rse6','rse7','rse8','rse9','rse10','rse11','rse12','rse14','0024_1B','0024_2B','0024_3B','0024_4'],0.38,0.4
;make histogram of snr for cluster members
;input maskarr : array of mask names
;      zmin    ; minimum z for members
;      zmax    ; maximum z for members

  goodfit  = []
  goodspec = []
  sn       = []
  snfit    = []
  mass     = []
  z        = []
  for i=0,n_elements(maskarr)-1 do begin
     mask = maskarr[i]
     directory = '/scr2/nichal/workspace2/sps_fit/data/'+mask+'/'
     sciencefits = directory+'sps_fit.fits.gz'
     scienceall = mrdfits(sciencefits,1,/silent)
     
     good = where(scienceall.z gt zmin and scienceall.z lt zmax, cgood)
     goodfit  = [goodfit,scienceall(good).goodfit]
     goodspec = [goodspec,scienceall(good).good]
     sn       = [sn,scienceall(good).sn]
     snfit    = [snfit,scienceall(good).snfit]
     mass     = [mass,scienceall(good).logmstar]
     z        = [mass,scienceall(good).z]
  endfor

  wsn = where(snfit gt 0, cwsn)
  if cwsn gt 0 then sn(wsn) = snfit(wsn)


  set_plot,'ps'
  !p.multi = [0,2,2]
  !p.font = 0
  psname='snr_hist.eps'
  device, filename = psname,xsize = 20,ysize = 16, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  binsnr = 2
  plothist,sn,xhist,yhist,bin=binsnr,xtitle='SN',ytitle='number of galaxies'
  goodfitpersn = fltarr(n_elements(xhist))
  for i=0,n_elements(xhist)-1 do begin
     insn = where(sn ge binsnr*i and sn lt binsnr*(i+1),cinsn)
     if cinsn gt 0 then goodfitpersn[i] = total(goodfit(insn))
  endfor
  cgplot,xhist,goodfitpersn,color=fsc_color('red'),psym=14,/overplot
  al_legend,['all galaxies','good fit'],psym=[15,14],color=fsc_color(['black','red']),/right,/top
  xyouts,0.5,0.9,'total galaxies: '+strcompress(string(n_elements(goodfit), format='(I)'),/rem),/normal
  xyouts,0.5,0.8,'total good spectra: '+strcompress(string(total(goodfit), format='(I)'),/rem),/normal
  xyouts,0.5,0.7,'masks used: ',/normal
  for i=0,ceil(n_Elements(maskarr)/7.)-1 do begin
     xyouts,0.5,0.65-i*0.05,strjoin(maskarr[i*7:i*7+6<n_Elements(maskarr)-1],', '),/normal
  endfor
  plothist,mass,mxhist,yxhist,/autobin,xstyle=4,ystyle=4,/nodata
  plothist,mass,mxhist,yxhist,bin=0.5,xrange=[8,12],xtitle='log(mass)',ytitle='number of galaxies'
  wgoodfit = where(goodfit eq 1)
  plothist,mass(wgoodfit),bin=0.5,color=fsc_color('red'),/overplot
  device,/close
  stop
end
