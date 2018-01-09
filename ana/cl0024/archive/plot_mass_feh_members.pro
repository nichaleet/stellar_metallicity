pro plot_mass_feh_members,mask,clustername


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
psname=clustername+'_FeH_mass.eps'
device, filename = psname,xsize = 15,ysize = 10, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
ctotgals=0
ctotbadgals=0

for nm=0,nmasks-1 do begin
   sciencefits = directory+mask[nm]+'/sps_fit.fits.gz'
   scienceall = mrdfits(sciencefits,1,/silent)
   
   goodobj = where(scienceall.goodfit eq 1 and scienceall.zfit ge 0.38 and scienceall.zfit le 4.0 and scienceall.good eq 1, cgals)
   badobj  = where(scienceall.goodfit eq 0 and scienceall.zfit ge 0.38 and scienceall.zfit le 4.0 and scienceall.feh ne -999. and scienceall.good eq 1, cbadgals)
   mass = scienceall[goodobj].logmstar
   feh  = scienceall[goodobj].feh
   feherr = scienceall[goodobj].feherr
   zfit   = scienceall[goodobj].zfit
   objnames = scienceall[goodobj].objname
   
   badmass = scienceall[badobj].logmstar
   badfeh  = scienceall[badobj].feh
   badfeherr = scienceall[badobj].feherr
   badzfit   = scienceall[badobj].zfit
   badobjnames = scienceall[badobj].objname

   zcat = [0.1,0.2,0.3,0.4,0.55,0.7,1.0]
   color  = ['deep pink','orange red','gold','dark green','blue','purple']
   if nm eq 0 then plot,mass,feh,xtitle='Log(M/Msun)',ytitle='[Fe/H]',/nodata,xrange=[8.,11.5],xstyle=1,yrange=[-1.,0.3],ystyle=1
;   count=1
   for i=2,2 do begin
      sel = where(zfit gt zcat(i) and zfit le zcat(i+1),csel)
      if csel gt 0 then begin
         print,mask[nm]+' at z:'+strtrim(string(zcat[i],format='(F3.1)'),2)
         oploterror,mass(sel),feh(sel),feherr(sel),color=fsc_color(color(i)),psym=1,errcolor=fsc_color(color(i)),thick=2
         cgplot,mass(sel),feh(sel),color=fsc_color(color(i)),psym=14,/overplot
         ctotgals = ctotgals+csel
      endif

      sel = where(badzfit gt zcat(i) and badzfit le zcat(i+1),csel)
      if csel gt 0 then begin
         oploterror,badmass(sel),badfeh(sel),badfeherr(sel),color=fsc_color('gray'),psym=1,errcolor=fsc_color('gray'),thick=2
         cgplot,badmass(sel),badfeh(sel),color=fsc_color('gray'),psym=14,/overplot
         ctotbadgals = ctotbadgals+csel
      endif
   endfor
endfor

;Add Choi data
;(mass,FeH,FeH_err)
z01 = {zlow:0.1,zhigh:0.2,mass:[9.9,10.2,10.4,10.7,11.0],Feh:[-0.05,-0.06,-0.01,-0.03,0.02],Feherr:[0.04,0.02,0.01,0.01,0.01]}
z02 = {zlow:0.2,zhigh:0.3,mass:[10.2,10.5,10.7,11.0,11.3],Feh:[-0.08,-0.06,-0.03,-0.01,-0.05],Feherr:[0.04,0.02,0.01,0.01,0.02]}
z03 = {zlow:0.3,zhigh:0.4,mass:[10.5,10.8,11.0,11.3],Feh:[-0.11,-0.05,-0.02,-0.03],Feherr:[0.03,0.01,0.01,0.02]}
z04 = {zlow:0.3,zhigh:0.4,mass:[10.8,11.1,11.3],Feh:[-0.07,-0.04,-0.05],Feherr:[0.02,0.01,0.02]}
z06 = {zlow:0.3,zhigh:0.4,mass:[10.9,11.0,11.3],Feh:[-0.15,-0.02,-0.05],Feherr:[0.07,0.03,0.03]}
oploterror,z01.mass,z01.feh,z01.feherr,color=fsc_color(color(0)),linethick=2,errcolor=fsc_color(color(0))
cgplot,z01.mass,z01.feh,psym=2,color=fsc_color(color(0)),/overplot
oploterror,z02.mass,z02.feh,z02.feherr,color=fsc_color(color(1)),linethick=2,errcolor=fsc_color(color(1))
cgplot,z02.mass,z02.feh,psym=2,color=fsc_color(color(1)),/overplot
oploterror,z03.mass,z03.feh,z03.feherr,color=fsc_color(color(2)),linethick=2,errcolor=fsc_color(color(2))
cgplot,z03.mass,z03.feh,psym=2,color=fsc_color(color(2)),/overplot
oploterror,z04.mass,z04.feh,z04.feherr,color=fsc_color(color(3)),linethick=2,errcolor=fsc_color(color(3))
cgplot,z04.mass,z04.feh,psym=2,color=fsc_color(color(3)),/overplot
oploterror,z06.mass,z06.feh,z06.feherr,color=fsc_color(color(4)),linethick=2,errcolor=fsc_color(color(4))
cgplot,z06.mass,z06.feh,psym=2,color=fsc_color(color(4)),/overplot

;Labelling
zarr_str = strarr(n_elements(zcat)-1)
for nz=0,n_elements(zcat)-2 do zarr_Str[nz]=strtrim(string(zcat[nz],format='(F3.1)'),2)+'<z<'+strtrim(string(zcat[nz+1],format='(F3.1)'),2)
al_Legend,[zarr_str,'Cl 0024 data','Choi et al. (2014)'],psym=[15,15,15,15,15,15,14,2],color=[color,'black','black'],box=0,thick=2,charsize=1,symsize=[1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.],/right,/bottom,font=0
USERSYM, [-1,1], [0, 0]
oplot,[11.45],[-0.88],psym=2,thick=2
oplot,[11.4],[-0.88],psym=8
;al_Legend,['Choi et al. (2014)'],psym=-2,position=[10.6,-1.95]
device,/close
print,'total of '+strtrim(string(ctotgals),2)+' galaxies were plotted'
print,'total of '+strtrim(string(ctotbadgals),2)+' bad galaxies were plotted'

stop

end
