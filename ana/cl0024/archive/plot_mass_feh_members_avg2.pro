pro plot_mass_feh_members_avg2,mask
;plot_mass_feh_members_avg2,['rse1','rse2','rse3','rse4','rse5','rse6','rse7','rse8','rse9','rse10','rse11','rse12','rse14','0024_1B','0024_2B','0024_3B','0024_4']

directory = '/scr2/nichal/workspace2/sps_fit/data/'
if ~keyword_set(mask) then begin
   maskname= strsplit(file_search(directory+'00*',/test_directory),'/',/extract)
   nmasks  = n_elements(maskname)
   nsub    = n_Elements(maskname[0])
   mask    = strarr(nmasks)
   for i=0,nmasks-1 do mask[i] = (maskname[i])[nsub-1]
endif

nmasks = n_elements(mask)
set_plot,'ps'
!p.multi = [0,1,1]
!p.font = 0
psname='allcl0024_FeH_mass_avg2.eps'
device, filename = psname,xsize = 15,ysize = 10, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
ctotgals=0
mass = []
feh  = []
feherr = []
zfit = []
objnames = []
goodfit = []
for nm=0,nmasks-1 do begin
   sciencefits = directory+mask[nm]+'/sps_fit.fits.gz'
   scienceall = mrdfits(sciencefits,1,/silent)
   
   goodobj = where(scienceall.good eq 1 and scienceall.zfit ge 0.38 and scienceall.zfit le 0.4, cgals)
   mass = [mass,scienceall[goodobj].logmstar]
   feh  = [feh,scienceall[goodobj].feh]
   feherr = [feherr,scienceall[goodobj].feherr]
   zfit   = [zfit,scienceall[goodobj].zfit]
   objnames = [objnames,scienceall[goodobj].objname]
   goodfit = [goodfit,scienceall[goodobj].goodfit]
endfor
plot,mass,feh,xtitle='Log(M/Msun)',ytitle='[Fe/H]',/nodata,xrange=[9.5,11.5],xstyle=5,yrange=[-.8,0.3],ystyle=5
zcat = [0.1,0.2,0.3,0.4,0.55,0.7,1.0]
color  = ['deep pink','gold','orange red','dark green','blue','purple']
zi = 2
;   count=1
massarr = [9.5,9.75,10.,10.25,10.55,11.]
nmass = n_elements(massarr)-1

sel = where(goodfit eq 1,csel)
marr = fltarr(nmass)+1./0
marr_lo = fltarr(nmass)+1./0
marr_hi = fltarr(nmass)+1./0
feharr  = fltarr(nmass)+1./0
fehstdarr = fltarr(nmass)+1./0
feherrarr = fltarr(nmass)+1./0
if csel gt 0 then begin
   zmass = mass(sel)
   zfeh  = feh(sel)
   zfeherr = feherr(sel)
   for mi=0,nmass-1 do begin
      marr(mi) = (massarr[mi]+massarr[mi+1])*0.5
      if mi eq nmass-1 then msel = where(zmass gt massarr[mi] and zfeh gt -1.,cmsel) else msel = where(zmass gt massarr[mi] and zmass lt massarr[mi+1] and zfeh gt -1.,cmsel)
      print, cmsel
      if cmsel gt 1 then begin
         meanerr,zfeh(msel),zfeherr(msel),fehmean,sigmam,sigmad,sigmas
         feharr(mi) = fehmean
         feherrarr(mi) = sigmam
         fehstdarr(mi) = sigmas
         meanerr,zmass(msel),zfeherr(msel),mmean,a,b,c
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
endif
   ;Add Choi data
;(mass,FeH,FeH_err)
vsym,5,/star
z01 = {zlow:0.1,zhigh:0.2,mass:[9.9,10.2,10.4,10.7,11.0],Feh:[-0.05,-0.06,-0.01,-0.03,0.02],Feherr:[0.04,0.02,0.01,0.01,0.01]}
z02 = {zlow:0.2,zhigh:0.3,mass:[10.2,10.5,10.7,11.0,11.3],Feh:[-0.08,-0.06,-0.03,-0.01,-0.05],Feherr:[0.04,0.02,0.01,0.01,0.02]}
z03 = {zlow:0.3,zhigh:0.4,mass:[10.5,10.8,11.0,11.3],Feh:[-0.11,-0.05,-0.02,-0.03],Feherr:[0.03,0.01,0.01,0.02]}
z04 = {zlow:0.3,zhigh:0.4,mass:[10.8,11.1,11.3],Feh:[-0.07,-0.04,-0.05],Feherr:[0.02,0.01,0.02]}
z06 = {zlow:0.3,zhigh:0.4,mass:[10.9,11.0,11.3],Feh:[-0.15,-0.02,-0.05],Feherr:[0.07,0.03,0.03]}

oploterror,z01.mass,z01.feh,z01.feherr,color=fsc_color(color(0)),linethick=2,errcolor=fsc_color(color(0))
oplot,z01.mass,z01.feh,psym=8,color=fsc_color(color(0)),symsize=1.5
oploterror,z02.mass,z02.feh,z02.feherr,color=fsc_color(color(1)),linethick=1.5,errcolor=fsc_color(color(1))
oplot,z02.mass,z02.feh,psym=8,color=fsc_color(color(1)),symsize=1.5
oploterror,z04.mass,z04.feh,z04.feherr,color=fsc_color(color(3)),linethick=1.5,errcolor=fsc_color(color(3))
oplot,z04.mass,z04.feh,psym=8,color=fsc_color(color(3)),symsize=1.5
oploterror,z06.mass,z06.feh,z06.feherr,color=fsc_color(color(4)),linethick=1.5,errcolor=fsc_color(color(4))
oplot,z06.mass,z06.feh,psym=8,color=fsc_color(color(4)),symsize=1.5

oploterror,z03.mass,z03.feh,z03.feherr,color=fsc_color(color(2)),linethick=1.5,errcolor=fsc_color(color(2))
oplot,z03.mass,z03.feh,psym=8,color=fsc_color(color(2)),symsize=1.5

   ;plot individual measurements
sel_lax = indgen(N_elements(mass))
sel = sel_lax
oploterror,mass(sel),feh(sel),feherr(sel),color=fsc_color('gray'),psym=1,errcolor=fsc_color('gray'),thick=2
cgplot,mass(sel),feh(sel),color=fsc_color('gray'),psym=14,/overplot

      ;plot average
cgerrplot,marr,feharr-feherrarr,feharr+feherrarr,color=fsc_color(color(zi))
cgerrplot,feharr,marr_lo,marr_hi,color=fsc_color(color(zi)),/horizontal
cgplot,marr,feharr,color=(color(zi)),psym=14,symsize=1.5,/overplot
ctotgals = ctotgals+csel

;Labelling
zarr_str = strarr(n_elements(zcat)-1)
for nz=0,n_elements(zcat)-2 do zarr_Str[nz]=strtrim(string(zcat[nz],format='(F3.1)'),2)+'<z<'+strtrim(string(zcat[nz+1],format='(F3.1)'),2)
vsym,5,/star
al_Legend,[zarr_str[0:4],'Cl 0024 data','Choi et al. (2014)'],psym=[15,15,15,15,15,14,8],color=[color[0:4],'black','black'],box=0,thick=2,charsize=1,symsize=[1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5],/right,/bottom,font=0
plot,mass,feh,xtitle='Log(M!D*!N/M!Isun!N)',ytitle='[Fe/H]',/nodata,xrange=[9.5,11.5],xstyle=1,yrange=[-1.,0.3],ystyle=1,/noerase

device,/close
print,'total of '+strtrim(string(ctotgals),2)+' galaxies were plotted'
stop

end
