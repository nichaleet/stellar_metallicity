pro linearfunc,x,a,f,pder
       f=a[0]+a[1]*x
       pder = [[replicate(1.,n_elements(x))],[x]]
end

function linfit_fixedslope_mc,x,y,dyarr,probdyarr,pars,sigma,chisq=chisq
;fit linear function with a fixed slope via monte carlo method. useful when the probability distribution of y is 
;nog gaussian
;INPUT
;x and y are arrays of what you fit (arrays of ngals size)
;yarr and probyarr are ngrid*nobjs grid (each yarr[*,i] and proyxarr is the probability for ith object)
   nmc = 1000. ;do 1000 timesi
   guess_param = pars ;fit parameters

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
      dy(i) = 0.5*(interpol(curyarr,cumprob[*,i],0.84)-interpol(curyarr,cumprob[*,i],0.16))
   endfor

   ;run over nmc loops to get parameters
   for i=0,nmc-1 do begin   
      ysamp = fltarr(nobjs)
      randomarr = randomu(seed,nobjs)
      for j=0,nobjs-1 do begin
         ysamp(j) = y(j)+interpol(dyarr[*,j],cumprob[*,j],randomarr(j))
      endfor
      pnow = guess_param
      bestfit = curvefit(x,ysamp,1./dy^2,pnow,function_name='linearfunc',status=status,fita=[1,0])
      if status ne 0 then stop,'curvefit failed'
      paramarr[*,i] = pnow
   endfor

   pars = mean(paramarr,dimension=2)
   sigma = stddev(paramarr,dimension=2)
   sigma[1] = 0     
   yfit = pars[0]-pars[1]*x
   if arg_present(chisq) then begin
      chisq = total((yfit-y)^2/dy^2)
   endif
   return,yfit
end
