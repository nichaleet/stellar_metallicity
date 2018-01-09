pro exam_sps
  set_plot,'x'
  f = file_search('/scr2/nichal/workspace2/fsps-3.0/SSP/SSP_Padova_MILES_Kroupa_*.spec',count=cfile)
  len = strpos(f,'.out')-strpos(f,'zmet')-4.
  filenum = intarr(cfile)
  for i=0,cfile-1 do filenum[i] = fix(strmid(f[i],strpos(f[i],'zmet')+4,len[i]))
  f = f(sort(filenum))
  spsall = sps_read_spec(f)
  nsize = size(spsall,/dimensions)
  nage  = nsize(0)
  nmet  = nsize(1)
;  for j=16,20 do for i=0,10 do begin
;  for j=18,18 do for i=0,10 do begin
  for j=10,30 do for i=3,3 do begin
;  for j=10,10 do for i=0,0 do begin
     sps = spsall[i*2+70,j]
     lambda = sps.lambda
     spsspec = sps.spec

     w = where(lambda gt 3700 and lambda lt 5500, c)
     if c lt 25 then message, 'Not enough pixels.'
     lambda = lambda[w]
     spsspec = spsspec[w]
     
     readcol, '/scr2/nichal/workspace2/sps_fit/lines.txt', linestart, lineend, linetype, format='D,D,A,X', /silent, comment='#'
     contmask = bytarr(n_elements(lambda))+1
     for ii=0,n_elements(linestart)-1 do begin
        w = where(lambda ge linestart[ii] and lambda le lineend[ii], c)
        if c gt 0 then contmask[w] = 0
     endfor
     won = where(contmask eq 1, con)
     if con lt 25 then message, 'Not enough pixels.'
     bkpt = slatec_splinefit(lambda[won], spsspec[won], coeff, bkspace=150, upper=3, lower=3, /silent,/everyn)
     if bkpt[0] eq -1 then message, 'Could not fit a spline to spsspec.'
     cont = slatec_bvalu(lambda, bkpt, coeff)
     spsspec /= cont
     plot,lambda,spsspec,title=strtrim(string(sps.zmet),2)+' '+strtrim(string(sps.agegyr),2)
     stop
  endfor

  ;plot fix age, change metallicity
  set_plot,'ps'
  psname = 'exam_sps_fixage.eps'
  device, filename = psname,xsize = 20,ysize = 18, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  !p.multi = [0,1,1]
  !p.font=0
  !p.charsize=0
  spssel=reform(spsall[80,16:20])
  color=['purple','blue','green','orange','red']
  zarr = fltarr(5)
  for i=0,n_elements(spssel)-1 do begin
     sps = spssel[i]
     lambda = sps.lambda
     spsspec = sps.spec

     w = where(lambda gt 3700 and lambda lt 5500, c)
     if c lt 25 then message, 'Not enough pixels.'
     lambda = lambda[w]
     spsspec = spsspec[w]
     if i eq 0 then plot,lambda,spsspec/median(spsspec),xtitle='lambda',/nodata
     oplot,lambda,spsspec/median(spsspec),color=fsc_color(color[i])
     zarr[i] = sps.zmet
  endfor
  al_legend,'[Fe/H]='+string(zarr,format='(F4.1)'),psym=15,color=color
  xyouts,3600,1.5,'age(gyr) = '+string(sps.agegyr,format='(F3.1)')
  device,/close

  ;plot fix metallicity, change age
  set_plot,'ps'
  psname = 'exam_sps_fixzmet.eps'
  device, filename = psname,xsize = 20,ysize = 18, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
  !p.multi = [0,1,1]
  !p.font=0
  !p.charsize=0
  selage = findgen(5)*4+70
  spssel=reform(spsall[selage,18])
  color=['purple','blue','green','orange','red']
  agearr = fltarr(5)
  for j=0,n_elements(spssel)-1 do begin
     i = 4-j
     sps = spssel[i]
     lambda = sps.lambda
     spsspec = sps.spec

     w = where(lambda gt 3700 and lambda lt 5500, c)
     if c lt 25 then message, 'Not enough pixels.'
     lambda = lambda[w]
     spsspec = spsspec[w]
     med = median(spsspec)
     ;med = 1.4e-15
     if i eq 4 then plot,lambda,spsspec/med,xtitle='lambda',ytitle='flux/median(flux)',/nodata
     oplot,lambda,spsspec/med,color=fsc_color(color[i])
     agearr[i] = sps.agegyr
  endfor
  al_legend,'Age(gyr)='+string(agearr,format='(F4.1)'),psym=15,color=color
  xyouts,3600,1.5,'[Fe/H]= '+string(sps.zmet,format='(F4.1)')
  device,/close
  stop

end
