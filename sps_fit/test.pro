pro test ;compare sps spectra old version and new version
  clight = 299792.458
  vdisp  =500. 
  f = file_search('/scr2/nichal/workspace2/fsps/SSP/*.spec')
  sps = sps_read_spec(f)
  f = file_search('/scr2/nichal/workspace2/fsps/SSP_old/*.spec')
  sps_old = sps_read_spec(f)
  
  dim = size(sps_old.spec,/dimensions)
  !p.multi=[0,1,2]
  setplot,5
  r = 'g'
  i = 0
  while i le dim(1)-1 and r ne 'q' do begin
     for j=0,dim(2)-1 do begin
        zmet = sps_old[i,j].zmet
        age = sps_old[i,j].agegyr
        loc = where(sps.zmet eq zmet and sps.agegyr eq age,c)
        if c ne 0 then begin
           lambdaold= sps_old[i,j].lambda
           specold  = sps_old[i,j].spec
           w = where(lambdaold gt 0 and lambdaold lt 10000,c)
           lambdaold=lambdaold[w]
           specold  =specold[w] 
           bkpt = slatec_splinefit(lambdaold, specold, coeff, bkspace=330, upper=5, lower=1, /silent,/everyn)
           if bkpt[0] eq -1 then message, 'Could not fit a spline to spsspec.'
           contold = slatec_bvalu(lambdaold, bkpt, coeff)
           specold_div = specold/contold
           specold_div = smooth_gauss_wrapper(lambdaold, specold_div, lambdaold, vdisp/clight/2.35*lambdaold)
           
           lambda= sps[loc].lambda
           spec  = sps[loc].spec
           w = where(lambda gt 0 and lambda lt 10000,c)
           lambda=lambda[w]
           spec  =spec[w]
           bkpt = slatec_splinefit(lambda, spec, coeff, bkspace=330, upper=5, lower=1, /silent,/everyn)
           if bkpt[0] eq -1 then message, 'Could not fit a spline to spsspec.'
           cont = slatec_bvalu(lambda, bkpt, coeff)
           spec_div = spec/cont
           spec_div = smooth_gauss_wrapper(lambda, spec_div, lambda, vdisp/clight/2.35*lambda)
           title='z:'+strtrim(string(zmet),2)+' age:'+strtrim(string(age),2)
           plot,lambdaold,specold_div,title=title
           oplot,lambda,spec_div,color=55
           plot,lambdaold,specold
           oplot,lambda,spec,color=55
           oplot, lambdaold,contold, linestyle=2,color=100
           oplot,lambda,cont,linestyle=2,color=100
           r=get_kbrd()
        endif
     endfor
     i = i+1
  endwhile
end
