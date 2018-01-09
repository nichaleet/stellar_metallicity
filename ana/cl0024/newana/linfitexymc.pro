function linfitexymc,x,y,dx,dy,dxarr,dyarr,probdxdyarr,sigma,chisq=chisq,yfit=yfit,plot=plot
;fit linear via monte carlo method. useful when the probability distribution of y is 
;not gaussian
;INPUT
;x and y are arrays of what you fit (arrays of ngals size)
;yarr and probyarr are ngrid*nobjs grid (each yarr[*,i] and proyxarr is the probability for ith object)
   nmc = 1000. ;do 1000 timesi
   paramarr = fltarr(2,nmc)

   dimen = size(probdxdyarr,/dimensions)
   ngridx = dimen(0)
   ngridy = dimen(1)
   nobjs = dimen(2)
   ;getting cumulative probability
   cumprob = dblarr(ngridx,ngridy,nobjs)
   dxgrid = fltarr(ngridx,ngridy,nobjs)
   dygrid = fltarr(ngridx,ngridy,nobjs)
   for gg=0,nobjs-1 do begin
      volume=0.d
      for ii=0,ngridx-1 do begin
         for jj=0,ngridy-1 do begin
            if ii eq 0 then deltax=abs(dxarr[ii+1,gg]-dxarr[ii,gg]) else deltax = abs(dxarr[ii,gg]-dxarr[ii-1,gg])
            if jj eq 0 then deltay=abs(dyarr[jj+1,gg]-dyarr[jj,gg]) else deltay = abs(dyarr[jj,gg]-dyarr[jj-1,gg])
            cumprob[ii,jj,gg] = volume
            volume = volume+deltax*deltay*probdxdyarr[ii,jj,gg]
         endfor
      endfor
      if abs(volume-1.) gt 0.01 then stop,'oops'
      dxgrid[*,*,gg] = rebin(dxarr[*,gg],ngridx,ngridy)
      dygrid[*,*,gg] = transpose(rebin(dyarr[*,gg],ngridy,ngridx))
   endfor
   ;run over nmc loops to get parameters
   if keyword_set(plot) then begin
      set_plot,'x'
      !p.multi=[0,1,1]
   endif
   for mm=0,nmc-1 do begin   
      xsamp = fltarr(nobjs)
      ysamp = fltarr(nobjs)
      randomarr = randomu(seed,nobjs)
      for j=0,nobjs-1 do begin
         probdiff = cumprob[*,*,j]-randomarr(j)
         boobee = min(abs(probdiff),locmin)
         minprobdiff = probdiff(locmin)
         probdiff(locmin) = max(abs(probdiff))
         boobee = min(abs(probdiff),loc2ndmin)
         secminprobdiff = probdiff(loc2ndmin)
         ;xsamp(j) = x(j)+interpol([dxgrid(locmin),dxgrid(loc2ndmin)],[minprobdiff,secminprobdiff],0.)
         ;ysamp(j) = y(j)+interpol([dygrid(locmin),dygrid(loc2ndmin)],[minprobdiff,secminprobdiff],0.)
         xsamp(j) = x(j)+dxgrid(locmin)
         ysamp(j) = y(j)+dygrid(locmin)
         ;if randomu(seed) gt 0.99 then print,dxgrid(locmin),dygrid(locmin)
      endfor
      fitexy,xsamp,ysamp,A,B,x_sig=dx,y_sig=dy,sigma_a_b
      ;param = linfit(xsamp,ysamp)
      ;a = param(0)
      ;b = param(1)
      if keyword_set(plot) and mm mod 250 eq 0 then begin
         plot,xsamp,ysamp,psym=1
         oplot,!x.crange,A+B*!x.crange,color=255
         ;stop
      endif
      paramarr[*,mm] = [A,B]
   endfor

   pout = median(paramarr,dimension=2)
   sigma = stddev(paramarr,dimension=2)
   A=pout(0)
   B=pout(1)
   if keyword_set(plot) then begin
      set_plot,'x'
      !p.multi=[0,1,2]
      plothist,paramarr[0,*],xtitle='A (const)'
      oplot,[A,A],!y.crange,color=255
      oplot,[A-sigma[0],A-sigma[0]],!y.crange,color=255,linestyle=2
      oplot,[A+sigma[0],A+sigma[0]],!y.crange,color=255,linestyle=2
      plothist,paramarr[1,*],xtitle='B (slope)'
      oplot,[B,B],!y.crange,color=255
      oplot,[B-sigma[1],B-sigma[1]],!y.crange,color=255,linestyle=2
      oplot,[B+sigma[1],B+sigma[1]],!y.crange,color=255,linestyle=2

      print,pout,sigma
      stop,'check the plots'
   endif
        
   if arg_present(chisq) or arg_present(yfit) then begin
      yfit = pout[0]+pout[1]*x
      chisq = total((yfit-y)^2/dy^2)
   endif
   return,pout
end
