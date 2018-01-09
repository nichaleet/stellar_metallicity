pro plotspec_panel,clusterarr,objnumarr,lambdalim=lambdalim,ylim=ylim,nofit=nofit,smooth=smooth,annotearr=annotearr
;plotspec_panel,['0024_2B','Cl1604_D2','Cl1604_D1','Cl1604_B'],[56,17,15,2],lambdalim=[3800,5500],smooth=[0,0,1,0],nofit=[0,1,1,1],annotearr=['Cl0024 z=0.4 DEIMOS','Cl1604 z=0.9, Good DEIMOS & Good MOSFIRE','Cl1604 z=0.9, Bad DEIMOS & Good MOSFIRE','Cl1604 z=0.9, Good LRIS & Bad MOSFIRE']
nspec = n_elements(clusterarr)
if ~keyword_set(lambdalim) then lambdalim = minmax(science.lambda)
if ~keyword_Set(ylim) then ylim = [-.5,2.5]
if ~keyword_set(nofit) then nofit = [1,1,1,1] 
if ~keyword_set(smooth) then smooth = intarr(nspec)
ysize = 4*nspec+2
xleft = 0.1
xright = 0.95
panelsize = 0.85/nspec

psname = '/scr2/nichal/workspace2/sps_fit/plotspec/panel_spec.eps'
set_plot,'ps'
device, filename = psname,xsize = 18,ysize = ysize, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
!p.multi=[0,1,1]
!p.font = 0
for i=0,nspec-1 do begin
   cluster= clusterarr[i]
   objnum = objnumarr[i]
   if keyword_set(annotearr) then annote = annotearr[i]

   file= '/scr2/nichal/workspace2/sps_fit/data/'+cluster+'/sps_fit.fits.gz'
   scienceall=mrdfits(file,1)
   science = scienceall[objnum-1]
   
   t = round(-1*ts_diff(science.fitmask, 1))
   wstart = where(t eq 1, cstart)+1
   wend = where(t eq -1, cend)
   if science.fitmask[0] eq 1 then begin
     if cstart eq 0 then begin
        wstart = 0
     endif else begin
        wstart = [0, wstart]
     endelse
     cstart += 1
  endif
   if science.fitmask[n_elements(t)-1] eq 1 then begin
      if cend eq 0 then begin
         wend = n_elements(t)-1
      endif else begin
        wend = [wend, n_elements(t)-1]
     endelse
      cend += 1
   endif
   if cstart ne cend then message, 'There are a different number of starting and ending fitmask wavelengths.'
  
   if science.zfit ne 0 and finite(science.zfit) then z=science.zfit else z=science.zspec
   goodrange = where(science.lambda/(1d + z) gt lambdalim[0] and science.lambda/(1d + z) lt lambdalim[1])
   lambda = science.lambda(goodrange)/(1d + z)
   spec = science.contdiv(goodrange) ;/science.spscont(goodrange)
   ivar = science.contdivivar(goodrange)

   if smooth[i] eq 1 then begin
      mask = science.fitmask(goodrange)
      if i eq 2 then bad = where(spec gt 1.25 or spec lt 0.75 and lambda gt 5000. and lambda lt 5300.)
     ;bad = where(lambda gt 5000. and lambda lt 5500.)
      good = where(mask eq 1 and lambda gt 5000.)
      meandiv = stdev(spec(good))
      spec(Bad) = (randomn(seed,n_elements(bad))*meandiv)+1.
                                ;spec(bad) = smooth(spec(bad),3,/nan)
   endif
   yhigh  = 0.95-panelsize*i
   ylow   = 0.95-panelsize*(i+1)
   if i eq 0 then plot, lambda,spec, xrange=lambdalim, yrange=ylim, xstyle=5, ystyle=5, color=fsc_color('black'),position=[xleft,ylow,xright,yhigh],/nodata
   if i gt 0 then plot, lambda,spec, xrange=lambdalim, yrange=ylim, xstyle=5, ystyle=5, color=fsc_color('black'),position=[xleft,ylow,xright,yhigh],/nodata,/noerase

   if cstart eq 0 or cend eq 0 then message, 'There are no fitmask regions.', /info else begin
      for j=0,cstart-1 do begin
         x = ([science.lambda[wstart[j]], science.lambda[wstart[j]], science.lambda[wend[j]], science.lambda[wend[j]]]/(1d + science.zspec) > lambdalim[0]) < lambdalim[1]
         y = [ylim[0], ylim[1], ylim[1], ylim[0]]
         polyfill, x, y, color=fsc_color('light cyan'),/data
      endfor
   endelse
   oplot, [-1d6, 1d6], [0.0, 0.0], color=fsc_color('pink')
   oplot, [-1d6, 1d6], [1.0, 1.0], color=fsc_color('pale green')
   oplot, lambda,spec, color=fsc_color('black') 
   
   readcol,'/scr2/nichal/workspace2/telluric/telluric.mask', tellstart, tellend, format='D,D', /silent, comment='#'
   tellthick = [5, 2, 5, 2, 2,2]
   linewaves = [2798.0, 3646.00, 3727.425, 3750.15, 3770.63, 3797.90, 3835.39, 3868.71, 3888.65, 3889.05, 3933.663, 3967.41, 3968.468, 3970.07, 4101.76, 4305.05, 4340.47, 4861.33, 4958.92, 5006.84, 5167.321, 5172.684, 5183.604, 5875.67, 5889.951, 5895.924, 6300.30, 6548.03, 6562.80, 6583.41, 6678.152, 6716.47, 6730.85]
   linenames = ['MgII', 'Hbreak', '[OII]', 'H12', 'H11', 'H10', 'H9', ' ', ' ', 'H8', 'CaH', ' ', 'CaK', ' ', 'Hd', 'CH', 'Hg', 'Hb', '[OIII]', '[OIII]', ' ', 'Mgb', ' ', 'HeI', 'NaD', 'NaD', '[OI]', '[NII]', 'Ha', '[NII]', 'HeI', '[SII]', '[SII]']
   linecolors = ['blue', 'black', 'blue', 'black', 'black', 'black', 'black', 'blue', 'blue', 'black', 'red', 'blue', 'red', 'black', 'black', 'red', 'black', 'black', 'blue', 'blue', 'red', 'red', 'red', 'blue', 'red', 'red', 'blue', 'blue', 'black', 'blue', 'blue', 'blue', 'blue']
   
   n = n_elements(tellstart)
   for j=0,n-1 do begin
      oplot, [(tellstart)[j], (tellend)[j]] / (1d + science.zspec), 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(tellthick)[j]
   endfor
   n = n_elements(linewaves)
   for j=0,n-1 do begin
      if (linewaves)[j] le !X.CRANGE[0] or (linewaves)[j] ge !X.CRANGE[1] then continue
      oplot, [(linewaves)[j], (linewaves)[j]], [0.06*!Y.CRANGE[0]+0.94*!Y.CRANGE[1], 0.02*!Y.CRANGE[0]+0.98*!Y.CRANGE[1]], color=fsc_color((linecolors)[j])
      xyouts, (linewaves)[j]+0.002*(!X.CRANGE[1]-!X.CRANGE[0]), 0.07*!Y.CRANGE[0]+0.93*!Y.CRANGE[1], (linenames)[j], orientation=90, alignment=1, color=fsc_color((linecolors)[j])
   endfor
   
   if i lt nspec-1 then plot, lambda,spec, xrange=lambdalim, yrange=ylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), /nodata, /noerase,position=[xleft,ylow,xright,yhigh],xtickformat='(A1)',yticks=6,ytickv=[0,0.5,1.,1.5,2,2.5]
   if i eq nspec-1 then plot, lambda,spec, xrange=lambdalim, yrange=ylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!3rest wavelength (!sA!r!u !9o!n)!3', /nodata, /noerase,position=[xleft,ylow,xright,yhigh],yticks=6,ytickv=[0,0.5,1.,1.5,2,2.5]
   
  if science.feh gt -10 and science.age gt 0.0 and total(science.spsspec gt 0.0) and nofit[i] eq 0 then oplot, science.lambda/(1d + science.zspec), science.spsspec, color=fsc_color('red')
  
  if keyword_set(annotearr) then xyouts,3850,-0.3,annote
endfor

xyouts,0.05*!d.x_size,0.45*!d.y_size,'!3flux (normalized)!3',orientation=90.,/device,charsize=1.

device,/close
stop
end
