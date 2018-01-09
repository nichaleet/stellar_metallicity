pro meanmc,x,dxarr,probdxarr,xmean,sigmam,sigmad,sigmas,plot=plot  
;find the mean via monte carlo method. useful when the probability distribution of x is not gaussian
;x is what you want to find mean (nobjs)
;dxarr and probdxarr are ngrid*nobjs grid (each dxarr[*,i] and prodbxarr is the probability for ith object)
   nmc = 1000. ;do 1000 times
   xmean_arr = fltarr(nmc)
   xsigma_arr= fltarr(nmc)  

   dimen = size(dxarr,/dimensions)
   ngrid = dimen(0)
   nobjs = dimen(1)

   ;getting cumulative probability
   cumprob = fltarr(ngrid,nobjs)
   dx = fltarr(nobjs)
   for i=0,nobjs-1 do begin
      curdxarr = dxarr[*,i]
      curprob = probdxarr[*,i]
      for j=1,ngrid-1 do cumprob[j,i] = int_tabulated(curdxarr[0:j],curprob[0:j])  
   endfor 

   ;run over nmc loops to get mean
   for i=0,nmc-1 do begin
      xsamp = fltarr(nobjs)
      randomarr = randomu(seed,nobjs)
      for j=0,nobjs-1 do begin
         xsamp(j) = x(j)+interpol(dxarr[*,j],cumprob[*,j],randomarr(j)) 
      endfor
      xmean_arr(i) = mean(xsamp)
      xsigma_arr(i) = stdev(xsamp)
   endfor
   sigmad = stdev(x)
   xmean = mean(xmean_arr)
   sigmam = stdev(xmean_arr)
   sigmas = mean(xsigma_arr)
   if keyword_set(plot) then begin 
     set_plot,'x'
     plothist,xmean_arr,xtitle='x value'
     oplot,[xmean,xmean],!y.crange,color=255
     print,'mean = ', mean(x)
     print,'weighted mean =',xmean,'+/-',sigmam
     print,'stdev(weight,unweight) =',sigmas,sigmad
     stop,'not much, pause for plot, .cont to cont'
   endif
end
