pro checkchisq
  common sps_spec, sps, spsz, spsage
  if (size(sps))[1] eq 0 then spsspec = sps_interp(0.0, 5.0)
  
;read science
  dir = '/scr2/nichal/workspace2/sps_fit/data/'
  scifits = dir+'cluster1/sps_fit.fits.gz'
  clustername = 'cluster1'
  scienceall = mrdfits(scifits,1,/silent)
;make output directory
  outputdir = '/scr2/nichal/workspace2/SDSSana/chisq/'+clustername
;select only members
  member = where(scienceall.zfit gt 0.08 and scienceall.zfit lt 0.12,cmember)
  scienceall = scienceall(member)
;sort them by mass
  sort_pos = sort(scienceall.logmstar)
  scienceall = scienceall(sort_pos)
;start doing the chisq check
  n_obj = n_elements(scienceall)
  for i=0,n_obj-1 do begin
     science      = scienceall[i]
     feh    = science.feh
     feherr = science.feherr
     zfit   = science.zfit
     age    = science.age
     ageerr = science.ageerr
     vdisp  = science.vdisp
     vdisperr = science.vdisperr
     lmass    = science.logmstar

     ;below is the same as sps_fit::fit
     restlambda = science.lambda / (1d + science.zspec)
     znow = science.zspec
     reallambda = science.lambda
     nlambda = n_elements(reallambda)
     if min(science.dlam) gt 0.2 and max(science.dlam,/nan) lt 10.0 then begin
        dlam_all = science.dlam
     endif else begin
        specresfile = dir+'cluster1/specres_poly.sav'
        if file_test(specresfile) then begin
           restore, specresfile
           dlam_all = poly(science.lambda/1000 - 7.8, specres_poly) / 2.35
        endif else dlam_all = replicate(3.9/2.35, nlambda)
     endelse

    won = where(science.fitmask eq 1 and finite(science.contdiv) and finite(science.contdivivar) and science.contdivivar gt 0, con)
    
    xmp  = reallambda[won]
    ymp  = science.contdiv[won]
    dymp = (science.contdivivar[won])^(-0.5)
    dlam = dlam_all[won]

    ;make grids around the found_feh
    feh_arr = (findgen(50.)*0.02)-0.5+feh
    chisq_feh = fltarr(50)
    for j=0,49 do begin
       fehnow = feh_arr[j]
       if fehnow lt max(spsz) and fehnow gt min(spsz) then begin
          a = [fehnow,age,vdisp,zfit]
          curspec = get_sps_obs(xmp,a,dlam)
          chisq = total((curspec-ymp)^2/dymp^2)
          chisq_feh[j] = chisq/float(n_elements(curspec)-4.) 
       endif else chisq_feh[j]=1./0.
    endfor
   ;make grids around the found age
    age_arr = findgen(15)+0.1
    chisq_age = fltarr(15)

    for j=0,14 do begin
       agenow = age_arr[j]
       if agenow lt max(spsage) and agenow gt min(spsage) then begin
          a = [feh,agenow,vdisp,zfit]
          curspec = get_sps_obs(xmp,a,dlam)
          chisq = total((curspec-ymp)^2/dymp^2)
          chisq_age[j] = chisq/float(n_elements(curspec)-4.) 
       endif else chisq_age[j]=1./0.
    endfor

    ;check if it's a part of the bad clump
    badclump = (feh lt .538*lmass-6.311 and lmass lt 10.8) or (feh lt -0.6 and lmass gt 10.8)
    ;plot
    psname = outputdir+'/'+strtrim(string(i),2)+'_'+science.objname+'.eps'
    set_plot,'ps'
    !p.multi = [0,2,1]
    device, filename = psname,xsize = 15,ysize = 6, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
    ;plot feh
    goodchi = where(finite(chisq_feh))
    plot,feh_arr,chisq_feh,psym=1,/nodata,yrange=[min(chisq_feh(goodchi))-0.1,max(chisq_feh(goodchi))+0.1],xtitle='[Fe/H]',ytitle='Reduced Chi-square'
    minchi = !y.crange[0]
    maxchi = !y.crange[1]
    cgcolorfill,[feh-feherr,feh+feherr,feh+feherr,feh-feherr],[maxchi,maxchi,minchi,minchi],color=fsc_color('rose')
    ;oplot with histogram
    goodrange= where(scienceall.feh lt max(feh_arr) and scienceall.feh gt min(feh_arr))
    bin=0.1
    plothist,scienceall(goodrange).feh,xhist,yhist,/noplot,peak=1,bin=bin
    yhist = yhist*(maxchi-minchi)/2.+minchi
    FOR k =0 ,N_Elements(xhist)-1 do cgcolorfill, [xhist[k]-bin/2., xhist[k]-bin/2.,xhist[k]+bin/2., xhist[k]+bin/2.], [minchi,yhist[k], yhist[k],minchi],color=fsc_color('cyan')
    vline,feh
    oplot,feh_arr,chisq_feh,psym=1
    ;make mark for badclump
    if badclump then al_legend,[''],psym=15,colors=['red'],/right
    ;plot age
    goodchi = where(finite(chisq_age))
    plot,age_arr,chisq_age,psym=1,/nodata,yrange=[min(chisq_age(goodchi))-0.1,max(chisq_age(goodchi))+0.1],xtitle='Age(Gyr)'
    minchi = !y.crange[0]
    maxchi = !y.crange[1]
    cgcolorfill,[age-ageerr,age+ageerr,age+ageerr,age-ageerr],[maxchi,maxchi,minchi,minchi],color=fsc_color('rose')
   ;oplot with histogram
    goodrange= where(scienceall.age lt max(age_arr) and scienceall.age gt min(age_arr))
    bin=1.
    plothist,scienceall(goodrange).age,xhist,yhist,/noplot,peak=1,bin=bin
    yhist = yhist*(maxchi-minchi)/2.+minchi
    FOR k =0 ,N_Elements(xhist)-1 do cgcolorfill, [xhist[k]-bin/2., xhist[k]-bin/2.,xhist[k]+bin/2., xhist[k]+bin/2.], [minchi,yhist[k], yhist[k],minchi],color=fsc_color('cyan')
    
    vline,age
    oplot,age_arr,chisq_age,psym=1
                                ;make mark for badclump
    if badclump then al_legend,[''],psym=15,colors=['red'],/right
    xyouts,0.7,0.9,science.logmstar,/device
    device,/close
    print, 'Finished doing', i
 endfor

end
