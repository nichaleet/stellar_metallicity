pro mzr,x,par,f,pder
  m0 = par[0]
  gamma = par[1]
  z0 = par[2]
  mass = x
  f = z0-alog10(1.+10.^(m0-x)*gamma)
  pder = fltarr(n_elements(x),n_elements(par)) ; no value returned.
end

pro plot_mass_feh_zform,mask,clustername
;e.g. plot_mass_feh_zform,['rse1','rse2','rse3','rse4','rse5','rse6','rse7','rse8','rse9','rse10','rse11','rse12','rse14','0024_1B','0024_2B','0024_3B','0024_4'],'cl0024'
 directory = '/scr2/nichal/workspace2/sps_fit/data/'
 if ~keyword_set(mask) then begin
    maskname= strsplit(file_search(directory+'00*',/test_directory),'/',/extract)
    nmasks  = n_elements(maskname)
    nsub    = n_Elements(maskname[0])
    mask    = strarr(nmasks)
    for i=0,nmasks-1 do mask[i] = (maskname[i])[nsub-1]
 endif

 nmasks = n_elements(mask)
;read data
 mass = []
 feh  = []
 feherr = []
 zfit = []
 age  = []
 objnames  = []

 for nm=0,nmasks-1 do begin
    sciencefits = directory+mask[nm]+'/sps_fit.fits.gz'
    scienceall = mrdfits(sciencefits,1,/silent)
   
    goodobj= where(scienceall.good eq 1 and scienceall.goodfit eq 1 and scienceall.zfit ge 0.38 and scienceall.zfit le 4.0, cgals)
    mass   = [mass,scienceall[goodobj].logmstar]
    feh    = [feh,scienceall[goodobj].feh]
    feherr = [feherr,scienceall[goodobj].feherr]
    zfit   = [zfit,scienceall[goodobj].zfit]
    age    = [age,scienceall[goodobj].age]
    objnames = [objnames,scienceall[goodobj].objname]
 endfor
    ageform  = (galage(zfit,1000.)/1.e9-age)>0
;plot
 ctotgals = 0
 set_plot,'ps'
 !p.multi = [0,1,1]
 !p.font = 0
 psname=clustername+'_FeH_mass_ageform.eps'
 device, filename = psname,xsize = 15,ysize = 10, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color

 agecat = [0,3.223,5.747,7.164,14] ;age of univ when gal formed correspond to z = [infinity,3,2,1,0.7,0]
 zcat = [1000,2,1,0.7,0.4]

 color  = ['darkred','orange','dark green','blue','black']
 ploterror,mass,feh,feherr,psym=1,errcolor=fsc_Color('lightgray'),xtitle='Log(M/Msun)',ytitle='[Fe/H]',xrange=[9.,11.5],xstyle=1,yrange=[-0.8,0.2]
 vsym,5,/fill
 nage = n_elements(Agecat)-1
 statusarr = intarr(nage)
 m0arr = dblarr(nage,2)
 z0arr = dblarr(nage,2)
 gammaarr = dblarr(nage,2)
 chisqarr = dblarr(nage)
 xplot = findgen(40)/20.+9.
 for i=0,nage-1 do begin
    sel = where(ageform gt agecat(i) and ageform le agecat(i+1) and mass gt 0,csel)
    ;if i eq nage-1 then sel = where(ageform gt agecat(i) and ageform le agecat(i+1) and mass gt 0 and feh gt -0.4,csel)
    if csel gt 0 then begin
       oplot,mass(sel),feh(sel),psym=8,color=fsc_color(color(i))
       print,' at age:'+strtrim(string(agecat[i],format='(F3.1)'),2)+'Gyr -',csel
       ctotgals = ctotgals+csel
       a=[9.,0.55,-0.05]
       fita=[1,1,1]
       if i eq 0 then fita=[1,1,0] 
       fit = curvefit(mass(sel),feh(sel),feherr(sel),a,fita=fita,sigma,chisq=chisq,function_name='mzr',/noderivative,status=status)
       m0arr[i,*]=[a[0],sigma[0]]
       z0arr[i,*]=[a[2],sigma[2]]
       gammaarr[i,*]=[a[1],sigma[1]]
       statusarr(i)=status
       chisqarr(i) = chisq
       if status eq 0 then begin
          yplot = a[2]-alog10(1.+10.^(a[0]-xplot)*a[1])
          oplot,xplot,yplot,color=fsc_color(color(i))
       endif
       if status ne 0 then stop
    endif
 endfor

;Labelling
 zarr_str = strarr(n_elements(zcat)-1)
 for nz=0,n_elements(zcat)-2 do zarr_Str[nz]=strtrim(string(zcat[nz],format='(F3.1)'),2)+'<z$\tex_{form}$<'+strtrim(string(zcat[nz+1],format='(F3.1)'),2)
 zarr_str(0) = 'z$\tex_{form}$>'+strtrim(string(zcat[1],format='(F3.1)'),2)
 al_Legend,zarr_str,psym=15,color=color,box=0,thick=2,charsize=1,symsize=1.5,/right,/bottom,font=0
 device,/close
 
 print,'total of '+strtrim(string(ctotgals),2)+' galaxies were plotted'
 stop

end
