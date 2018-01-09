pro mzrcurve, x, a, f, pder
        f=a[0]-alog10(1.+10.^((x-a[1])*(-1.*a[2])))-8.9
        ; FEH = 12+log(O/H)-8.9 = z0-log[1+(M*/M0)^-gmma]-8.9
        ;pder=[[replicate(1.0, N_ELEMENTS(X))],[a[2]/(alog(10)*(a[2]*(a[1]-x)))],[(x-a[1])/(alog(10)*(a[2]*(a[1]-x)))]]
        pder=[[replicate(1.0, N_ELEMENTS(X))],[replicate(1.0, N_ELEMENTS(X))],[replicate(1.0, N_ELEMENTS(X))]]
end

function fitzahidmc,mass,feh,dfeharr,probdfeharr,sigma,chisq=chisq,guess_param=guess_param
;fit log function in zahid13 via monte carlo method. useful when the probability distribution of feh is 
;nog gaussian
;INPUT
;mass and feh are arrays of what you fit (arrays of ngals size)
;feharr and probfeharr are ngrid*nobjs grid (each feharr[*,i] and profehxarr is the probability for ith object)
   nmc = 1000. ;do 1000 timesi
   if ~keyword_set(guess_param) then guess_param = [8.8,9.7,0.6] ;fit parameters

   paramarr = fltarr(3,nmc)

   dimen = size(dfeharr,/dimensions)
   ngrid = dimen(0)
   nobjs = dimen(1)

   ;getting cumulative probability
   cumprob = fltarr(ngrid,nobjs)
   dfeh = fltarr(nobjs)
   for i=0,nobjs-1 do begin
      curdfeharr = dfeharr[*,i]
      curprob = probdfeharr[*,i]
      for j=1,ngrid-1 do cumprob[j,i] = int_tabulated(curdfeharr[0:j],curprob[0:j])
      dfeh(i) = 0.5*(interpol(dfeharr[*,i],cumprob[*,i],0.84)-interpol(dfeharr[*,i],cumprob[*,i],0.16))
   endfor

   ;run over nmc loops to get parameters
   for i=0,nmc-1 do begin   
      fehsamp = fltarr(nobjs)
      randomarr = randomu(seed,nobjs)
      for j=0,nobjs-1 do begin
         fehsamp(j) = feh(j)+interpol(dfeharr[*,j],cumprob[*,j],randomarr(j))
      endfor
      pnow = guess_param
      bestfit = curvefit(mass,fehsamp,1./dfeh^2,pnow,function_name='mzrcurve',status=status)
      if status ne 0 then stop,'curvefit failed'
      paramarr[*,i] = pnow
   endfor

   pout = mean(paramarr,dimension=2)
   sigma = stddev(paramarr,dimension=2)
        
   if arg_present(chisq) then begin
      yfit = pout[0]-alog10(1.+10.^((mass-pout[1])*(-1.*pout[2])))-8.9
      chisq = total((yfit-feh)^2/dfeh^2)
   endif
   return,pout
end
