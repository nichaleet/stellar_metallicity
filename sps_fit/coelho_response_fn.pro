function coelho_response_fn,feh,age,afearr,vdisp,redshift
  common coelho_spec, c07, c07feh, c07age, c07afe
  common get_sps, dlam, dataivar, datalam, wonfit, contmask, normalize,rest    
  common get_coelho,coelho_normalize

  clight = 299792.458
;create response function with different [alpha/Fe] values at fixed feh, age
  base_struct = coelho_interp(feh,age,0)
  lambda = base_struct.lambda
  base_spec = base_struct.spec
  base_spec = smooth_gauss_wrapper(lambda*(redshift+1.), base_spec, datalam, sqrt(dlam^2+(vdisp/clight*datalam)^2))
  if coelho_normalize eq 1 then begin
     won = where(contmask eq 1, con)
     if con lt 25 then message, 'Not enough pixels.'
     innvar = dataivar[won]/(median(base_spec[won]))^2
     bkpt = slatec_splinefit(lambda[won]/(1.+redshift), base_spec[won], coeff, invvar=invvar, bkspace=100, upper=3, lower=3, /silent) 
     if bkpt[0] eq -1 then message, 'Could not fit a spline to spsspec.'
     cont = slatec_bvalu(lambda/(1.+redshift), bkpt, coeff)
     base_spec /= cont
  endif

  nafe = n_elements(afearr)
  for i=0,nafe-1 do begin
     spec_struct = coelho_interp(feh,age,afearr[i])
     lambda = spec_struct.lambda
     spec = spec_struct.spec
     spec = smooth_gauss_wrapper(lambda*(redshift+1.), spec, datalam, sqrt(dlam^2+(vdisp/clight*datalam)^2))

     if coelho_normalize eq 1 then begin
        innvar = dataivar[won]/(median(spec[won]))^2
        bkpt = slatec_splinefit(lambda[won]/(1.+redshift), spec[won], coeff, invvar=invvar, bkspace=100, upper=3, lower=3, /silent) 
        if bkpt[0] eq -1 then message, 'Could not fit a spline to spsspec.'
        cont = slatec_bvalu(lambda/(1.+redshift), bkpt, coeff)
        spec /= cont
     endif

     dspec = spec/base_spec

     if i eq 0 then begin 
        astr = {feh:feh,age:age,afe:afearr[i],vdisp:vdisp,redshift:redshift,dlam:dlam,lambda:datalam,dspec:dspec}
        str = replicate(astr,nafe)
     endif else begin
        str[i].afe = afearr[i]
        str[i].dspec = dspec
     endelse
  endfor

return, str

end
