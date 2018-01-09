pro plot_feh_alpha,mask
;plot_feh_alpha,['rse1','rse2','rse3','rse4','rse5','rse6','rse7','rse8','rse9','rse10','rse11','rse12','rse14','0024_1B','0024_2B','0024_3B','0024_4']
;plot_feh_alpha,['all_cl0024']
directory = '/scr2/nichal/workspace2/sps_fit/data/'
if ~keyword_set(mask) then begin
   maskname= strsplit(file_search(directory+'00*',/test_directory),'/',/extract)
   nmasks  = n_elements(maskname)
   nsub    = n_Elements(maskname[0])
   mask    = strarr(nmasks)
   for i=0,nmasks-1 do mask[i] = (maskname[i])[nsub-1]
endif

nmasks = n_elements(mask)
for i=0,1 do begin
   mass = []
   feh  = []
   feherr = []
   zfit = []
   objnames = []
   goodfit = []
   afe = []
   afeupper = []
   afelower = []
   goodsky = []
   age = []
   ageerr = []
   mass = []
   for nm=0,nmasks-1 do begin
	sciencefits = directory+mask[nm]+'/sps_fit.fits.gz'
	scienceall = mrdfits(sciencefits,1,/silent)
	goodobj = where(scienceall.goodfit+scienceall.good eq 2 and scienceall.zfit ge 0.38 and scienceall.zfit le 4.0, cgals)
        badobj  = where(scienceall.goodfit+scienceall.good eq 1 and scienceall.zfit ge 0.38 and scienceall.zfit le 4.0 and scienceall.feh ne -999., cbadgals)
        if i eq 0 then obj = goodobj
        if i eq 1 then obj = badobj
	mass = [mass,scienceall[obj].logmstar]
	feh  = [feh,scienceall[obj].feh]
	feherr = [feherr,scienceall[obj].feherr]
	zfit   = [zfit,scienceall[obj].zfit]
	objnames = [objnames,scienceall[obj].objname]
	goodfit = [goodfit,scienceall[obj].goodfit]
	afe = [afe,scienceall[obj].alphafe]
	afeupper = [afeupper,scienceall[obj].alphafeupper]
	afelower = [afelower,scienceall[obj].alphafelower]
	goodsky = [goodsky,scienceall[obj].goodsky]
	age = [age,scienceall[obj].age]
	ageerr = [ageerr,scienceall[obj].ageerr]
        mass = [mass,scienceall[obj].logmstar]
   endfor
   str = {mass:mass,feh:feh,feherr:feherr,zfit:zfit,objnames:objnames,goodfit:goodfit,afe:afe,afeupper:afeupper,afelower:afelower,goodsky:goodsky,age:age,ageerr:ageerr}
   if i eq 0 then goodstr = str
   if i eq 1 then badstr  = str 
endfor

set_plot,'ps'
!p.multi = [0,1,2]
!p.font = 0
psname='feh_alpha.eps'
device, filename = psname,xsize = 15,ysize = 15, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
icolor = ['darkgreen','gray']
for i=1,0,-1 do begin
    if i eq 0 then str=goodstr else str=badstr
    normal = where(str.afe lt 0.4 and str.afe gt 0.,cnormal)
    upperlim = where(str.afe le 0.,cupperlim)
    lowerlim = where(str.afe ge 0.4,clowerlim)
    if i eq 1 then plot,str.mass,str.afe,xtitle='Log(M/Msun)',ytitle='[Alpha/Fe]',xrange=[9.5,11.5],/nodata,yrange=[-0.1,0.5]
    if cnormal gt 0 then begin
	nolowerlim = where(str.afelower(normal) lt -10,cnolowerlim)
	if cnolowerlim gt 0 then str.afelower(normal(nolowerlim)) = 2.*str.afe(normal(nolowerlim))-str.afeupper(normal(nolowerlim))
        noupperlim = where(str.afeupper(normal) gt 10,cnoupperlim)
        if cnoupperlim gt 0 then str.afeupper(normal(noupperlim)) = 2.*str.afe(normal(noupperlim))-str.afelower(normal(noupperlim))
	cgerrplot,str.mass(normal),str.afelower(normal),str.afeupper(normal),color=icolor[i],psym=1
    endif
    goodsky = where(str.goodsky eq 1,cgoodsky,complement=badsky,ncomplement=cbadsky)
    vsym,4,/fill,rot=45
    if cgoodsky gt 0 then oplot,str.mass(goodsky),str.afe(goodsky),psym=8,color=fsc_color(icolor[i])
    if cbadsky gt 0 then cgplot,str.mass(badsky),str.afe(badsky),psym=16,color=icolor[i],/overplot 
    plotsym,1,2,thick=2,color=fsc_color(icolor[i])	
    if cupperlim gt 0 then oplot,str.mass(upperlim),str.afe(upperlim),psym=8
    plotsym,2,2,thick=2,color=fsc_color(icolor[i])
    if clowerlim gt 0 then oplot,str.mass(lowerlim),str.afe(lowerlim),psym=8
endfor

for i=1,0,-1 do begin
    if i eq 0 then str=goodstr else str=badstr
    if i eq 1 then plot,str.mass,str.age,xtitle='Log(M/Msun)',ytitle='Age(Gyr)',xrange=[9.5,11.5],/nodata,yrange=[1,12],/ylog,ystyle=1
       oploterror,str.mass,str.age,str.ageerr,errcolor=fsc_color(icolor[i]),psym=1
    goodsky = where(str.goodsky eq 1,cgoodsky,complement=badsky,ncomplement=cbadsky)
    vsym,4,/fill,rot=45
    if cgoodsky gt 0 then oplot,str.mass(goodsky),str.age(goodsky),psym=8,color=fsc_color(icolor[i])
    if cbadsky gt 0 then cgplot,str.mass(badsky),str.age(badsky),psym=16,color=icolor[i],/overplot
endfor

;for i=0,1 do begin
;    if i eq 0 then str=goodstr else str=badstr
;    normal = where(str.afe lt 0.4 and str.afe ge 0.,cnormal)
;    upperlim = where(str.afe le 0.,cupperlim)
;    lowerlim = where(str.afe ge 0.4,clowerlim)
;    if i eq 0 then plot,str.feh,str.afe,xtitle='[Fe/H]',ytitle='[Mg/Fe]',/nodata,yrange=[-0.1,0.5],xrange=[-0.4,0.2]
;    if cnormal gt 0 then begin 
;    	cgerrplot,str.feh(normal),str.afelower(normal),str.afeupper(normal),color=icolor[i],psym=1
;    	cgerrplot,str.afe(normal),str.feh(normal)-str.feherr(normal),str.feh(normal)+str.feherr(normal),/horizontal,color=icolor[i]
;    endif
;    vsym,4,/fill,rot=45
;    oplot,str.feh,str.afe,psym=8,color=fsc_color(icolor[i])
;
;    if cupperlim gt 0 then begin
;        plotsym,1,2,thick=2,color=fsc_color(icolor[i])
;	oplot,str.feh(upperlim),str.afe(upperlim),psym=8
;	cgerrplot,str.afe(upperlim),str.feh(upperlim)-str.feherr(upperlim),str.feh(upperlim)+str.feherr(upperlim),/horizontal,color=icolor[i]
;    endif
;    if clowerlim gt 0 then begin
;        plotsym,2,2,thick=2,color=fsc_color(icolor[i])
;        oplot,str.feh(lowerlim),str.afe(lowerlim),psym=8
;        cgerrplot,str.afe(lowerlim),str.feh(lowerlim)-str.feherr(lowerlim),str.feh(lowerlim)+str.feherr(lowerlim),/horizontal,color=icolor[i]
;     endif
;endfor
device,/close

set_plot,'ps'
!p.multi = [0,1,1]
!p.font = 0
psname='alpha_age.eps'
device, filename = psname,xsize = 15,ysize = 10, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color

for i=1,0,-1 do begin
    if i eq 0 then str=goodstr else str=badstr
    normal = where(str.afe lt 0.4 and str.afe gt 0.,cnormal)
    upperlim = where(str.afe le 0.,cupperlim)
    lowerlim = where(str.afe ge 0.4,clowerlim)
    if i eq 1 then plot,str.age,str.afe,xtitle='Age(Gyr)',ytitle='[Alpha/Fe]',xrange=[1,11.5],/nodata,yrange=[-0.1,0.5]
    if cnormal gt 0 then begin
        nolowerlim = where(str.afelower(normal) lt -10,cnolowerlim)
        if cnolowerlim gt 0 then str.afelower(normal(nolowerlim)) = 2.*str.afe(normal(nolowerlim))-str.afeupper(normal(nolowerlim))
        noupperlim = where(str.afeupper(normal) gt 10,cnoupperlim)
        if cnoupperlim gt 0 then str.afeupper(normal(noupperlim)) = 2.*str.afe(normal(noupperlim))-str.afelower(normal(noupperlim))
        cgerrplot,str.age(normal),str.afelower(normal),str.afeupper(normal),color=icolor[i],psym=1
    endif
    cgerrplot,str.afe,str.age-str.ageerr,str.age+str.ageerr,/horizontal,color=icolor[i]
    goodsky = where(str.goodsky eq 1,cgoodsky,complement=badsky,ncomplement=cbadsky)
    vsym,4,/fill,rot=45
    if cgoodsky gt 0 then oplot,str.age(goodsky),str.afe(goodsky),psym=8,color=fsc_color(icolor[i])
    if cbadsky gt 0 then cgplot,str.age(badsky),str.afe(badsky),psym=16,color=icolor[i],/overplot
    plotsym,1,2,thick=2,color=fsc_color(icolor[i])
    if cupperlim gt 0 then oplot,str.age(upperlim),str.afe(upperlim),psym=8
    plotsym,2,2,thick=2,color=fsc_color(icolor[i])
    if clowerlim gt 0 then oplot,str.age(lowerlim),str.afe(lowerlim),psym=8
endfor
	
device,/close
stop
end
