pro sps_iterproc, funcname, p, iter, fnorm, functargs=functargs, parinfo=pi, quiet=quiet, dof=dof
  common sps_iterproc, contiter
  if iter gt 1 then print, contiter, p[0], p[1], p[2],p[3], fnorm/dof, dof, format='(I3,2X,D6.3,1X,D5.2,2X,D6.1,2x,D6.3,1X,D8.3,1X,I4)'
end

pro fitspschoi,science
  common sps_spec, sps, spsz, spsage
  common sps_iterproc, contiter
  common get_sps, dlam, dataivar, datalam, wonfit
  clight = 299792.458
  print,'Fitting ...'
  restlambda = science.lambda
  znow = science.zspec
  reallambda = science.lambda*(science.zspec+1.)
  nlambda = n_elements(reallambda)

  pi = replicate({value:0d, fixed:0, limited:[1,1], limits:[0.D,0.D], parname:'', mpprint:0, mpformat:'', step:0d, tied:''}, 4)
  pi.value = double([0.0, 5.0, 200.0,znow])
  pi[0].limits = minmax(spsz)
  pi[1].limits = minmax(spsage)
  pi[2].limits = [0.0, 1000.0]
  pi[3].limits = [-0.05,0.05]+znow
  pi.step = double([0.1, 0.5, 25.0,0.0002])
  pi.parname = ['    Z', '  age', 'vdisp','redshift']
  pi.mpformat = ['(D6.3)', '(D5.2)', '(D6.1)','(D6.3)']
  
  won = where(science.fitmask eq 1 and finite(science.contdiv) and finite(science.contdivivar) and science.contdivivar gt 0, con)
  if con lt 10 then begin
     pi.value = [-999d, -999d, -999d]
     perror = [-999d, -999d, -999d]
     science.spsspec = -999d
     goto, done
  endif
    
  xmp = reallambda[won]
  ymp = science.contdiv[won]
  dymp = (science.contdivivar[won])^(-0.5)
  wontofit = won
  zdiff = 1.0
  agediff = 1.0
  vdispdiff = 1.0
  redshfdiff = 1.0
  contiter = 0
  print, '* * * * * * * * * * * * * * * * * * * *'
  print, strtrim(science.objname, 2)
  print, '* * * * * * * * * * * * * * * * * * * *'
  print, '  i Z/Z_sun   age sigma_v  redhift    chi^2  DOF'
  print, '--- ------- ----- ------- ---------  -------- ----'

  while abs(zdiff) gt 0.001 or abs(agediff) gt 0.001 or abs(vdispdiff) gt 0.001 or abs(redshfdiff) gt 0.0001 do begin
     contiter++
     dlam = reallambda*350./clight/2.35
     dataivar = science.contdivivar
     datalam = reallambda
     wonfit = wontofit
     pars = mpfitfun('get_sps_obs_sdss', xmp, ymp, dymp, parinfo=pi, /nocatch, bestnorm=bestnorm, dof=dof, perror=perror, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=1000, status=status, yfit=ympfit, iterproc='sps_iterproc')
     zdiff = (pi[0].value-pars[0])/pi[0].value
     agediff = (pi[1].value-pars[1])/pi[1].value
     vdispdiff = (pi[2].value-pars[2])/pi[2].value
     redshfdiff = (pi[3].value-pars[3])/pi[2].value
     pi.value = pars
     ;stop
     dlam = restlambda*350./clight/2.35
     spsbestfit = get_sps_rest_sdss(restlambda, pi.value)
     bkpt = slatec_splinefit(restlambda[won], science.contdiv[won] / spsbestfit[won], coeff, invvar=science.contdivivar[won], bkspace=330, upper=3, lower=3, /silent)
     if bkpt[0] eq -1 then begin
        pi.value = [-999d, -999d, -999d]
        perror = [-999d, -999d, -999d]
        science.spsspec = -999d
        science.spscont = -999d ;nicha added this line
        break
     endif
     cont = slatec_bvalu(restlambda, bkpt, coeff)
     ymp = science.contdiv[won] / cont[won]
     dymp = (science.contdivivar[won])^(-0.5) / cont[won]
  endwhile
  pcor = covar*0.
  for i=0,3 do for j=0,3 do pcor[i,j]=covar[i,j]/sqrt(covar[i,i]*covar[j,j])
  science.pcor = pcor
  science.covar = covar
  science.spscont = cont        
  science.spsspec = spsbestfit
  print, ' '
    
  done:
  science.feh = pi[0].value
  science.feherr = perror[0]
  science.age = pi[1].value
  science.ageerr = perror[1]
  science.zfit = pi[3].value
  science.vdisp = pi[2].value
  science.vdisperr = perror[2]
 
end
