function linfunc, p
common linfunc_xy,x,y,dyarr,probdyarr,pin,fita
  if total(fita) ge 1 then begin
    param = pin
    param(where(fita ne 1)) = p
  endif else param = p
  
  model = param[0]+param[1]*x
  ndata = n_elements(x)
  probi = fltarr(ndata)
  for i=0,ndata-1 do begin
     curyarr = y(i)+dyarr[*,i]
     curprob = probdyarr[*,i]
     cur_delta_model = curyarr-model(i)
     dy = curyarr[1]-curyarr[0]
     probarr_i = gauss_smooth(curprob,param[2]/dy,/edge_truncate)
     probi(i) = interpol(probarr_i,cur_delta_model,0.)
    ; plot,cur_delta_model,curprob
    ; oplot,cur_delta_model,probarr_i,color=255
  endfor
  L = -1.*total(alog(probi))
  return,L
end

function intrinsic_scatter_linear_nongauss,xi,yi,p0in,dyarr_in,probdyarr_in,fita=fitain
common linfunc_xy,x,y,dyarr,probdyarr,pin,fita
;fita = 1 if fix the parameter
;p0 should be const, slope, intrinsic scatter
   x = xi
   y = yi
   dyarr = dyarr_in
   probdyarr = probdyarr_in
   pin = p0in
   if keyword_set(fitain) then begin
     fita = fitain
     wherefit = where(fita ne 1,cwherefit)
     p0 = pin(wherefit)
   endif else begin
      fita=lonarr(n_elements(p0))
      p0 = pin
   endelse
   pout= amoeba(1.0e-5,function_name='linfunc',scale=1.e-3,p0=p0,ncalls=ncalls)

   if total(fita) ge 1 then begin
      param = pin
      param(where(fita ne 1)) = pout
   endif else param=pout
   return,[param]
end
