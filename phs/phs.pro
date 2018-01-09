pro phs
files = file_search('*/sps_fit.fits.gz')
nclusters = n_elements(files)
clusternames=strarr(nclusters)
for i=0,nclusters-1 do  clusternames[i] = (strsplit(files[i],'/',/extract))[0]
color = ['purple','cyan','navy','darkgreen','green','yellow','orange','red','indianred','brown','gray'] 
path = '/path/to/your/cluster/folders/'
metal     = []
metal_err = []
mass      = []
mass_err = []
clusternumber = []

;2) Set the plot files
set_plot,'ps'
psname = 'mass_metallicity.eps'
device, filename = psname,xsize = 20,ysize = 15, $
      xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
!p.font=0
!p.charsize=1.5
;3) Begin the for loop to read data and plot them
firstplot = 1
for i=0,nclusters-1 do begin
    scinow = mrdfits(files[i],1)
    ngals   = n_Elements(scinow)
    if ngals gt 1 and i ne 4 then begin  
       badfeh = where(scinow.feh eq 0.198, cbadfeh)
       if cbadfeh gt 0 then scinow(badfeh).feherr = 0.1
       if cbadfeh gt 0 then scinow(badfeh).feh = 0.3+randomn(seed,cbadfeh)*0.1
       metal  = [metal,scinow.feh]
       metal_err = [metal_err,scinow.feherr]
       mass   = [mass,scinow.logmstar]
       clusternumber = [clusternumber,fltarr(ngals)+i]
       if firstplot eq 1 then begin
          ploterror,scinow.logmstar,scinow.feh,scinow.feherr,psym=1,errcolor=fsc_color(color[i]),xrange=[10.,12],yrange=[-0.5,0.4],/nodata
;read gallazzi
          glz = mrdfits('/scr2/nichal/workspace2/gallazzi_data.fits',1,/silent)
          inmass = where(glz.logmass gt min(!x.crange))
          inz    = where(glz.zmin gt min(!y.crange))
          firstpoint = interpol(glz.z,glz.logmass,[min(!x.crange)])
          lastpoint  = interpol(glz.logmass,glz.zmin,[min(!y.crange)])
          cgcolorfill,[min(!x.crange),min(!x.crange),glz.logmass(inmass),reverse(glz.logmass(inz)),lastpoint],[min(!y.crange),firstpoint,glz.zmax(inmass),reverse(glz.zmin(inz)),min(!y.crange)],color=fsc_Color('cornsilk')
          oplot,glz.logmass,glz.z,color=fsc_color('chartreuse')
          oplot,glz.logmass,glz.zmin,psym=0,linestyle=2,color=fsc_color('chartreuse')
          oplot,glz.logmass,glz.zmax,psym=0,linestyle=2,color=fsc_color('chartreuse')
          ploterror,scinow.logmstar,scinow.feh,scinow.feherr,psym=1,errcolor=fsc_color(color[i]),xtitle='log(M)',ytitle='[Fe/H]',xrange=[10.,12],yrange=[-0.5,0.4],/noerase
          firstplot = 0
       endif else oploterror,scinow.logmstar,scinow.feh,scinow.feherr,psym=1,errcolor=fsc_color(color[i])
    endif
endfor

;4) Find the average metallicity at each mass bin. For this we have to decide how many mass bins are there once we see the scatter plot.
 ;find massbin everage
  ;eachbin has 13 members sort them by mass
  morder = sort(mass)
  for i=0,7 do begin
     lastin = (i+1)*15-1 
     if lastin gt n_Elements(mass)-1 then lastin = n_elements(mass)-1
     in       = morder[i*15:lastin]
     lowmass = where(mass(in) lt 10,clowmass,complement=goodmass)
     if clowmass gt 0 then in = in(goodmass)
     massin   = mass(in)
     massave  = mean(massin)
     fehin    = metal(in)
     fehinerr = metal_err(in)
     fehave   = wmean(fehin,fehinerr,/nan)
     oploterror,[massave],[fehave],[stdev(massin)],[stdev(fehin)],errcolor=fsc_color('black'),errthick=2
     cgplot,[massave],[fehave],color=fsc_color('black'),/overplot,psym=14,symsize=2
     ;stop
  endfor
device,/close
stop
end
