pro snr_mass
dir  = '/scr2/nichal/workspace2/sps_fit/data/'
masks     = file_search(dir+'00*',/test_directory,/mark_directory)
nmasks    = n_elements(masks)
colors = ['purple','red','hot pink','dark green','blue', 'navy blue']
zarr = [0,0.3,0.6,1.]
namearr  = strarr(nmasks)
all_mass = []
all_snr  = []
all_z    = []
set_plot,'ps'
!p.multi = [0,2,2]
!p.font = 0
psname='SNR_mass.eps'
device, filename = psname,xsize = 20,ysize = 16, $
        xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
;Mass and SNR by mask
for nm=0,nmasks-1 do begin
   maskname = strsplit(masks[nm],'/',/extract)
   maskname = maskname[n_elements(maskname)-1]
   namearr[nm] = maskname
   sciencefile = file_search(masks[nm]+'sps_fit.fits.gz')
   science = mrdfits(sciencefile[0],1)
   goodmass = where(science.logmstar ne 0 and finite(science.logmstar),cgmass)
   if nm eq 0 then plot,science[goodmass].sn,science[goodmass].logmstar,psym=1,xtitle='SNR',ytitle='Log M*',/nodata,xrange=[5,25],yrange=[7,11.5],ystyle=1
   oplot,science[goodmass].sn,science[goodmass].logmstar,psym=1,color=fsc_color(colors[nm]),thick=2
   all_mass = [all_mass,science[goodmass].logmstar]
   all_snr   = [all_snr,science[goodmass].sn]
   all_Z   = [all_z,science[goodmass].zspec]
endfor
al_Legend,namearr,psym=1,color=colors,box=0,thick=2,charsize=1,symsize=1.5,/right,/bottom,font=0

;Mass and SNR by redshifts
colors=['navy blue','dark green','dark red']
zarr_str = strarr(3)
for nz = 0,n_elements(zarr)-2 do begin
   goodz = where(all_z gt zarr[nz] and all_z le zarr[nz+1])
   if nz eq 0 then plot,all_snr[goodz],all_mass[goodz],psym=1,xtitle='SNR',ytitle='Log M*',/nodata,xrange=[5,25],yrange=[7,11.5],ystyle=1
   oplot,all_snr[goodz],all_mass[goodz],psym=1,color=fsc_color(colors[nz]),thick=2
   zarr_str[nz] = strtrim(string(zarr[nz],format='(F3.1)'),2)+'<z<'+strtrim(string(zarr[nz+1],format='(F3.1)'),2)
endfor
al_Legend,zarr_str,psym=1,color=colors,box=0,thick=2,charsize=1,symsize=1.5,/right,/bottom,font=0

;Histogram of mass
goodsnr = where(all_snr gt 5. and all_mass gt 0)
plothist,all_mass[goodsnr],xtitle='Log M*',bin=0.5,ytitle='Number'
device,/close
stop
end
