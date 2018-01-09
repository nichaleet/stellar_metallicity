pro test_evol_doable
; to test if we can measure evolution with uncertainties of 0.1 with intrinsic scatter from z=0.4 to z=0
;setting
evol = 0.08 ;dex
slope = 0.2 ;dex per log mass
nsample = 70
intscat = 0.1;dex (intrinsic scatter)
uncert = 0.1 ;dex (measurement uncertainties)
himass = 11 ;solar mass
lomass = 9.5 
intercept1 = -2.2
intercept2 = intercept1-evol
nloop = 100

!p.multi=[0,1,1]
zscore=fltarr(nloop)
tscore=fltarr(nloop)
for i=0,nloop-1 do begin
   mass1 = randomu(seed,nsample)*(himass-lomass)+lomass
   mass2 = randomu(seed,nsample)*(himass-lomass)+lomass
   idealfeh1 = mass1*slope+intercept1+randomn(seed,nsample)*intscat
   idealfeh2 = mass2*slope+intercept2+randomn(seed,nsample)*intscat
   obsfeh1 = idealfeh1+randomn(seed,nsample)*uncert
   obsfeh2 = idealfeh2+randomn(seed,nsample)*uncert
   linpar1 = linfit(mass1,obsfeh1,sigma=sigma1,chisq=chi1,yfit=yfit1)
   linpar2 = linfit(mass2,obsfeh2,sigma=sigma2,chisq=chi2,yfit=yfit2)
   ;slope z score
   zscore(i) = (linpar1(1)-linpar2(1))/sqrt(sigma1(1)^2+sigma2(1)^2)
   ;intercept mean t test
   sse1 = total((yfit1-obsfeh1)^2)
   sse2 = total((yfit2-obsfeh2)^2)
   dof = nsample+nsample-4.
   s2_12 = (sse1+sse2)/dof ;pooled residual variance
   s_12 = sqrt(s2_12*(1./nsample+1./nsample));+(mean(mass1))^2/(variance(mass1)*(nsample-1.))+(mean(mass2))^2/(variance(mass2)*(nsample-1.))))
;   s_12 = sqrt(s2_12)
;   tscore(i) = (linpar1(0)-linpar2(0))/s_12

    tscore(i) = (mean(obsfeh1)-mean(obsfeh2))/sqrt(variance(obsfeh1)/nsample+variance(obsfeh2)/nsample)

   plot,mass1,obsfeh1,psym=4
   oplot,mass1,obsfeh1,psym=4,color=fsc_color('cyan')
   oplot,mass2,obsfeh2,psym=4,color=fsc_color('salmon')
   oplot,mass1,yfit1,psym=0,color=fsc_color('blue')
   oplot,mass2,yfit2,psym=0,color=fsc_color('red')
   ;wait,0.1
endfor

!p.multi=[0,1,2]
plothist,zscore,bin=0.1
plothist,tscore
oplot,[-2,-2],!y.crange,color=fsc_color('red')
oplot,[2,2],!y.crange,color=fsc_color('red')
help,where(abs(tscore) gt 2.65) ;3 sigma
stop
end
