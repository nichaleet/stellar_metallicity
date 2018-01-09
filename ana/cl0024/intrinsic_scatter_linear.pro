function linfunc, p
common linfunc_xy,x,y,dy
L = 0.5*total(alog((p[2])^2+dy^2))+0.5*total((y-p[0]-p[1]*x)^2/(dy^2+p[2]^2))
return,L
end

function intrinsic_scatter_linear,xi,yi,dyi,p0,mc=mc,yarr=yarr,probyarr=probyarr
;linear fit and find intrinsic scatter at the same time
;if select montecarlo (mc) then it'll also return uncertainties of the best fit params using montecarlo method
;minimize chisquare according to https://arxiv.org/pdf/physics/0511182.pdf equation 69
;inputs are xi(x values), yi(y values), dyi(measurement errors), p0 = best guest for [intercept,slope,intrinsic scatter)
common linfunc_xy,x,y,dy
x = xi
y = yi
dy = dyi

if keyword_set(mc) then begin
   nmc = 1000
   ndata = n_elements(x)
   param = fltarr(3,nmc)
   chisq = fltarr(nmc)

   if keyword_set(yarr) then begin
      ;getting cumulative probability
      dimen = size(yarr,/dimensions)
      ngrid = dimen(0)

      cumprob = fltarr(ngrid,ndata)
      dy = fltarr(ndata)
      for i=0,ndata-1 do begin
         curyarr = yarr[*,i]
         curprob = probyarr[*,i]
         for j=1,ngrid-1 do cumprob[j,i] = int_tabulated(curyarr[0:j],curprob[0:j])
         dy(i) = 0.5*abs(interpol(yarr[*,i],cumprob[*,i],0.84)-interpol(yarr[*,i],cumprob[*,i],0.16))
      endfor
   endif
   
   for i=0,nmc-1 do begin
      if keyword_set(yarr) then begin
         ynow = fltarr(ndata)
         randomarr = randomu(seed,ndata)
         for j=0,ndata-1 do begin
            ynow(j) = interpol(yarr[*,j],cumprob[*,j],randomarr(j))
         endfor
      endif else ynow = y+randomn(seed,ndata)*dyi
      y = ynow
      param[*,i] = amoeba(1.0e-5,function_name='linfunc',scale=1.e-3,p0=p0,ncalls=ncalls)
      chisq[i] = total((y-param[0,i]-param[1,i]*x)^2/(dy^2+param[2,i]^2))
   endfor

   pout1 = fltarr(3)
   pout2 = fltarr(3)
   for i=0,2 do begin
      pout1[i] = mean(param[i,*])
      pout2[i] = stdev(param[i,*])
   endfor
   return,[pout1,mean(chisq),pout2]

endif else begin
   param= amoeba(1.0e-5,function_name='linfunc',scale=1.e-3,p0=p0,ncalls=ncalls)
   chisq = total((y-param[0]-param[1]*x)^2/(dy^2+param[2]^2))
   return,[param,chisq]
endelse

end
