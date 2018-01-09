pro cannon,plot=plot,quadratic=quadratic
    filesps = 'cannon_sps_str.fits'
    ;cannon_makesps_struct,wlmin=3500,wlmax=7500,fileout=filesps
    sps = mrdfits(filesps,1)
    z = sps.zmet
    logage = alog10(sps.agegyr)
    nspec = n_elements(sps)
    if keyword_set(quadratic) then begin
       xarr = fltarr(4,nspec)
       xarr[0,*] = z
       xarr[1,*] = logage
       xarr[2,*] = z^2
       xarr[3,*] = logage^2
       nparam = 5
    endif else begin
       xarr = fltarr(2,nspec)
       xarr[0,*] = z
       xarr[1,*] = logage
       nparam = 3
    endelse
    spsz = z[uniq(z, sort(z))]
    spslogage = logage[uniq(logage, sort(logage))]  ;increment is 0.5 dex in log.. so should do it in log    
    lambda = sps[0].lambda
    cwl = n_elements(lambda)
    
    parr = fltarr(nparam,cwl)
    if keyword_set(plot) then set_plot,'x'
    for i=0,cwl-1 do begin
       flambda = sps.contdiv[i]
       p = regress(xarr,flambda,chisq=chisq,const=p0,correlation=corr,yfit=yfit,ftest=ftest,status=status)
       if keyword_set(plot) then begin
          !p.font = 0
          plot,yfit,flambda,psym=1,yrange=minmax(flambda),xrange=minmax(yfit),xtitle='yfit',$
               ytitle='f lambda'
          oplot,!x.crange,!x.crange,color=fsc_color('red')
          xyouts,0.1,0.8,'status: '+strtrim(string(status),2),/normal
          xyouts,0.1,0.75,'f value: '+sigfig(ftest,2),/normal
          xyouts,0.1,0.7,'param (z age): '+strjoin(sigfig(p,2),','),/normal
          aa=get_kbrd()
       endif
       parr[*,i]=[p0,reform(p)]
    endfor
    save,parr,filename='cannon_param.sav'
    set_plot,'ps'
    psname='cannon_param.eps'
    device, filename = psname,xsize = 30,ysize = 15, $
             xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
       !p.multi=[0,1,2]
       !p.font=0
       plot,lambda,parr[0,*],xtitle='lambda',ytitle='const term',xrange=[3500,7500],xstyle=1
       plot,lambda,parr[1,*],xtitle='lambda',ytitle='linear coeff term',/nodata,xrange=[3500,7500],xstyle=1
       oplot,lambda,parr[1,*],color=fsc_color('maroon')
       oplot,lambda,parr[2,*],color=fsc_color('navy')
       al_legend,['[Z/H]','log(age)'],psym=15,color=fsc_color(['maroon','navy']),box=0,/top,/right
       oplotlinenames
;       plot,lambda,parr[1,*]/parr[2,*],xtitle='lambda',ytitle='Z_coeff/Age_coeff',xrange=[3500,7500],xstyle=1,yrange=[-100,100]
;       oplotlinenames
    device,/close
stop
end
