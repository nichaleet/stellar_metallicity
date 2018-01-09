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
set_plot,'ps'
!p.multi = [0,1,1]
!p.font = 0
psname=clustername+'_FeH_mass_ageform.eps'
device, filename = psname,xsize = 15,ysize = 10, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
ctotgals=0
for nm=0,nmasks-1 do begin
   sciencefits = directory+mask[nm]+'/sps_fit.fits.gz'
   scienceall = mrdfits(sciencefits,1,/silent)
   
   goodobj= where(scienceall.good eq 1 and scienceall.goodfit eq 1 and scienceall.zfit ge 0.38 and scienceall.zfit le 4.0, cgals)
   mass   = scienceall[goodobj].logmstar
   feh    = scienceall[goodobj].feh
   feherr = scienceall[goodobj].feherr
   zfit   = scienceall[goodobj].zfit
   age    = scienceall[goodobj].age
   objnames = scienceall[goodobj].objname
   ageform  = galage(zfit,1000.)/1.e9-age
   ageform(where(ageform lt 0.)) = 0.

   agecat = [0,2.109,3.223,5.747,7.164,14] ;correspond to z = [infinity,3,2,1,0.7,0]
   zcat = [1000,3,2,1,0.7,0]

   color  = ['darkred','orange red','gold','dark green','blue','black']
   if nm eq 0 then ploterror,mass,feh,feherr,psym=1,xtitle='Log(M/Msun)',ytitle='[Fe/H]',/nodata,xrange=[9.,11.5],xstyle=1,title='SSP Fit'
;   count=1
   for i=0,n_elements(agecat)-2 do begin
      sel = where(ageform gt agecat(i) and ageform le agecat(i+1),csel)
      if csel gt 0 then begin
         print,mask[nm]+' at z:'+strtrim(string(agecat[i],format='(F3.1)'),2)
         oploterror,mass(sel),feh(sel),feherr(sel),color=fsc_color(color(i)),psym=1,errcolor=fsc_color(color(i)),thick=2
         ctotgals = ctotgals+csel
;         for j=0,csel-1 do begin
;            xyouts,mass[sel(j)]+0.05,feh[sel(j)],strtrim(string(count),2)
;            xyouts,[11.3],[-0.6-count*0.1],strtrim(string(count),2)+' '+objnames[sel(j)]
;            count=count+1
;         endfor
      endif
   endfor
endfor


;Labelling
zarr_str = strarr(n_elements(zcat)-1)
for nz=0,n_elements(zcat)-2 do zarr_Str[nz]=strtrim(string(zcat[nz],format='(I2)'),2)+'<z$\tex_{form}$>'+strtrim(string(zcat[nz+1],format='(I2)'),2)
zarr_str(0) = 'z$\tex_{form}$>'+strtrim(string(zcat[1],format='(I2)'),2)
al_Legend,zarr_str,psym=15,color=color,box=0,thick=2,charsize=1,symsize=1.5,/right,/bottom,font=0
device,/close

print,'total of '+strtrim(string(ctotgals),2)+' galaxies were plotted'
stop

end
