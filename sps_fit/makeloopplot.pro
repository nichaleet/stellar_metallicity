pro makeloopplot,lambda,data,fit,cont,nnan,nloop,id,chisq,pars,xmp,ymp
set_plot,'ps'
psname = 'examine/id_'+id+'_nloop'+strtrim(string(fix(nloop)),2)+'.eps'
device, filename = psname,xsize = 18,ysize = 9, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
!p.multi=[0,1,1]

xrange=[4500.,7000.]
plot,lambda,data,yrange=[0.5,1.3],xrange=xrange,xstyle=1, ystyle=1,color=fsc_color('black'),ytitle='!3flux!3',position=[0.1,0.7,0.9,0.96],xtickformat="(A1)",font=0
oplot,lambda[nnan],fit,color=fsc_color('red')

plot,xmp,ymp,xrange=xrange,yrange=[0.5,1.3],xstyle=1,ystyle=1,ytitle='ymp',position=[0.1,0.44,0.9,0.7],xtickformat='(A1)',/noerase,font=0
oplot,lambda[nnan],fit,color=fsc_color('red')

plot,lambda[nnan],data[nnan]/fit,xrange=xrange,yrange=[0.9,1.1],ytitle='data/fit',position=[0.1,0.18,0.9,0.44],xtitle='!3rest wavelength (!sA!r!u !9o!n)!3',xstyle=1,/noerase,ystyle=1,font=0
oplot,lambda,cont,color=fsc_color('blue')

mask = bytarr(n_elements(lambda))
mask[nnan] = 1
t = round(-1*ts_diff(mask, 1))
wstart = where(t eq 1, cstart)+1
wend = where(t eq -1, cend)
if mask[0] eq 1 then begin
   if cstart eq 0 then begin
      wstart = 0
   endif else begin
      wstart = [0, wstart]
   endelse
   cstart += 1
endif
if mask[n_elements(t)-1] eq 1 then begin
   if cend eq 0 then begin
      wend = n_elements(t)-1
   endif else begin
      wend = [wend, n_elements(t)-1]
   endelse
   cend += 1
endif
for i=0,cstart-1 do begin
   x = [lambda[wstart[i]],lambda[wend[i]]]
   y = !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0])+[0.,0.]
   oplot, x, y, color=fsc_color('lightcyan'), thick=5
endfor
xyouts,0.5,0.45,'[Fe/H]='+strtrim(string(pars[0],format='(F5.2)'),2)+', Age='+strtrim(string(pars[1],format='(F4.2)'),2)+', nloop ='+strtrim(string(fix(nloop)),2)+', !9'+string("143B)+'!x!u2!n='+strtrim(string(chisq,format='(F4.2)'),2),/normal,alignment=0.5,font=0
device,/close
set_plot,'x'

end
