pro mass_feh_zform_sdss
  
  ;read catalog
  readcol,'/scr2/nichal/workspace2/SDSSdata/gallazzi1/z_mstar_all.dat',plate,mjd,fiber,logmasslow,logmass,logmasshigh,FeHlow,feh,fehhigh,logagelow,logage,logagehigh,comment='#'
  age = 10.^(logage-9.)
  agelow = 10.^(logagelow-9.)
  agehigh = 10.^(logagehigh-9.)

  ageform = 12.3-age
  ageformlow = 12.3-agehigh 
  ageformhigh = 12.3-agelow
  ;median z is ~0.1 which is about 12.3 Gyr 

  agecat = [0,2.109,3.223,5.747,7.164,14] ;correspond to z = [infinity,3,2,1,0.7,0]
  zcat = [1000,3,2,1,0.7,0]

;  mcat=[9,10,10.25,10.4,10.6,10.7,10.8,10.9,11,11.1,11.2,11.3,11.5]
  mcat=[9,10,10.25,10.5,10.75,11,11.25,11.5]
  ;plot
  set_plot,'ps'
  !p.multi = [0,1,1]
  !p.font = 0
  psname='Gallazzi_FeH_mass_ageform.eps'

  device, filename = psname,xsize = 15,ysize = 10, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  color  = ['darkred','orange red','gold','dark green','blue','black']
  plot,logmass,feh,psym=1,xtitle='Log(M/Msun)',ytitle='[Fe/H]',/nodata,xrange=[9.,11.5],xstyle=1,yrange=[-1,0.5]
  ctotgals = 0

  nage = n_elements(agecat)-1
  nmass= n_Elements(mcat)-1
  feharr    = fltarr(nage,nmass)
  fehlowarr = fltarr(nage,nmass)
  fehhigharr= fltarr(nage,nmass)
  massarr    = fltarr(nage,nmass)
  masslowarr = fltarr(nage,nmass)
  masshigharr= fltarr(nage,nmass)
  vsym,5,/fill
  for i=0,nage-1 do begin
     zsel = where(ageform gt agecat(i) and ageform le agecat(i+1),czsel)
     print, 'z='+strtrim(string(zcat[i],format='(F3.1)'),2)+':'+strtrim(string(zcat[i+1],format='(F3.1)'),2)+' ngals =',czsel

     if czsel gt 0 then begin
        for j=0,nmass-1 do begin
           msel = where(logmass(zsel) gt mcat(j) and logmass(zsel) le mcat(j+1),csel)
           sel = zsel(msel)
           meanerr,logmass(sel),abs(logmass(sel)-logmasslow(sel)),meanmasslow,sigmamasslow,sigmad,errmasslow
           meanerr,logmass(sel),abs(logmass(sel)-logmasshigh(sel)),meanmasshigh,sigmamasshigh,sigmad,errmasshigh
           meanmass = wmean([meanmasslow,meanmasshigh],[sigmamasslow,sigmamasshigh])
           massarr(i,j) = meanmass
           masslowarr(i,j) = meanmass-errmasslow
           masshigharr(i,j) = meanmass+errmasshigh

           meanerr,feh(sel),abs(feh(sel)-fehlow(sel)),meanfehlow,sigmafehlow,sigmad,errfehlow
           meanerr,feh(sel),abs(feh(sel)-fehhigh(sel)),meanfehhigh,sigmafehhigh,sigmad,errfehhigh
           meanfeh = wmean([meanfehlow,meanfehhigh],[sigmafehlow,sigmafehhigh])
           feharr(i,j) = meanfeh
           fehlowarr(i,j) = meanfeh-errfehlow
           fehhigharr(i,j) = meanfeh+errfehhigh
           ctotgals = ctotgals+csel 
           if csel lt 10 then print,agecat(i),meanmass,csel
       endfor
     endif
  endfor

  ;plot
  cgerrplot,massarr,fehlowarr,fehhigharr,color='lightgray'
  cgerrplot,feharr,masslowarr,masshigharr,color='lightgray',/horizontal

  for i=0,nage-1 do begin
     oplot,massarr(i,*),feharr(i,*),color=fsc_color(color(i)),psym=8
  endfor

;Labelling
  zarr_str = strarr(n_elements(zcat)-1)
  for nz=0,n_elements(zcat)-2 do zarr_Str[nz]=strtrim(string(zcat[nz],format='(I2)'),2)+'>z$\tex_{form}$>'+strtrim(string(zcat[nz+1],format='(I2)'),2)
  zarr_str(0) = 'z$\tex_{form}$>'+strtrim(string(zcat[1],format='(I2)'),2)
  al_Legend,zarr_str,psym=15,color=color,box=0,thick=2,charsize=1,symsize=1.5,/left,/bottom,font=0
  device,/close

  print, ctotgals
  stop
end
