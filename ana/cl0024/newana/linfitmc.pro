function linfitmc,x,y,dyarr,probdyarr,sigma,chisq=chisq,yfit=yfit
;fit linear via monte carlo method. useful when the probability distribution of y is 
;not gaussian
;INPUT
;x and y are arrays of what you fit (arrays of ngals size)
;yarr and probyarr are ngrid*nobjs grid (each yarr[*,i] and proyxarr is the probability for ith object)
   nmc = 1000. ;do 1000 timesi

   paramarr = fltarr(2,nmc)

   dimen = size(dyarr,/dimensions)
   ngrid = dimen(0)
   nobjs = dimen(1)

   ;getting cumulative probability
   cumprob = fltarr(ngrid,nobjs)
   dy = fltarr(nobjs)
   for i=0,nobjs-1 do begin
      curyarr = dyarr[*,i]
      curprob = probdyarr[*,i]
      for j=1,ngrid-1 do cumprob[j,i] = int_tabulated(curyarr[0:j],curprob[0:j])
      dy(i) = 0.5*(interpol(dyarr[*,i],cumprob[*,i],0.84)-interpol(dyarr[*,i],cumprob[*,i],0.16))
   endfor

   ;run over nmc loops to get parameters
   for i=0,nmc-1 do begin   
      ysamp = fltarr(nobjs)
      randomarr = randomu(seed,nobjs)
      for j=0,nobjs-1 do begin
         ysamp(j) = y(j)+interpol(dyarr[*,j],cumprob[*,j],randomarr(j))
      endfor

      paramarr[*,i] = linfit(x,ysamp,measure_errors=dy)
   endfor

   pout = mean(paramarr,dimension=2)
   sigma = stddev(paramarr,dimension=2)
        
   if arg_present(chisq) or arg_present(yfit) then begin
      yfit = pout[0]+pout[1]*x
      chisq = total((yfit-y)^2/dy^2)
   endif
   return,pout
end
