; =================
pro sps_iterproc, funcname, p, iter, fnorm, functargs=functargs, parinfo=pi, quiet=quiet, dof=dof
    common sps_iterproc, contiter
    common toprint, agediff, zdiff
    if iter gt 1 then print, contiter, p[0], p[1], p[2],p[3], fnorm/dof, dof,abs(zdiff),abs(agediff),format='(I4,2X,D6.3,1X,D5.2,2X,D6.1,2x,D6.3,1X,D8.3,1X,I4,2X,D8.4,2X,D8.4)'
 end

pro get_vdisp, x, par, f, pder
common get_vdisp, dlamfit
  clight = 299792.458
  redshift= par[0]
  vdisp = par[1]
  Mgb1line = 5168.761*(1.+redshift)
  Mgb2line = 5174.125*(1.+redshift)
  Mgb3line = 5185.048*(1.+redshift)
  Mgb1dip = par[2]
  Mgb2dip = par[3]
  Mgb3dip = par[4]
  const = par[5]
  dlamMgb = median(dlamfit[where(x/(1.+redshift) gt 5160 and x/(1.+redshift) lt 5190.)])
  ww1 = (vdisp/clight/2.35*Mgb1line)^2+dlamMgb^2
  ww2 = (vdisp/clight/2.35*Mgb2line)^2+dlamMgb^2
  ww3 = (vdisp/clight/2.35*Mgb3line)^2+dlamMgb^2

  f = const -  Mgb1dip*exp(-0.5*(x-Mgb1line)^2/ww1)-Mgb2dip*exp(-0.5*(x-Mgb2line)^2/ww2)-Mgb3dip*exp(-0.5*(x-Mgb3line)^2/ww3)
  pder = fltarr(n_elements(x),n_elements(par)) ; no value returned.

end


pro sps_fit::indices, science, noredraw=noredraw, nostatusbar=nostatusbar
    nmc = 1000
    cah = 3933.663
    cahindex = dblarr(nmc+1) - 1d9
    gbandindex = dblarr(nmc+1) - 1d9

    restlambda = science.lambda / (1d + science.zspec)
    flux = science.contdiv
    ivar = science.contdivivar
    w = where(finite(flux) and finite(ivar) and ivar gt 0, n)
    if n lt 10 then return
    restlambda = restlambda[w]
    flux = flux[w]
    ivar = ivar[w]
    wcah = where(restlambda gt cah-9 and restlambda lt cah+9, ccah)
    wcahside = where((restlambda gt 3900 and restlambda lt 3916) or (restlambda gt 4015 and restlambda lt 4029), ccahside)
    wgband = where(restlambda gt 4285 and restlambda lt 4318, cgband)
    wgbandside = where(restlambda gt 4228 and restlambda lt 4273, cgbandside)
    for i=0,nmc do begin
        if i eq 0 then fluxi = flux else fluxi = flux + (ivar)^(-0.5)*randomn(seed, n)
        if ccah gt 10 and ccahside gt 10 then cahindex[i] = weightedmean(1.0 - fluxi[wcah], (ivar[wcah])^(-0.5)) - weightedmean(1.0 - fluxi[wcahside], (ivar[wcahside])^(-0.5))
        if cgband gt 10 and cgbandside gt 10 then gbandindex[i] = weightedmean(1.0 - fluxi[wgband], (ivar[wgband])^(-0.5)) - weightedmean(1.0 - fluxi[wgbandside], (ivar[wgbandside])^(-0.5))
    endfor
    science.cah = cahindex[0]
    w = where(cahindex[1:nmc] gt -1d8)
    science.caherr = stddev(cahindex[w+1])
    science.gband = gbandindex[0]
    w = where(gbandindex[1:nmc] gt -1d8)
    science.gbanderr = stddev(gbandindex[w+1])
end

pro sps_fit::vdisp_ppxf, science
  goodlam = where(science.lambda/(1.+science.zspec) gt 4800. and science.lambda/(1.+science.zspec) lt 5380., ngoodlam)
  
  ;;find dlam
  if min(science.dlam) gt 0.2 and max(science.dlam) lt 10.0 then begin
     dlam_all = science.dlam(goodlam)
  endif else begin
     specresfile = self.directory+'specres_poly.sav'
     if file_test(specresfile) then begin
        restore, specresfile
        dlam_all = poly(science.lambda/1000 - 7.8, specres_poly) / 2.35
     endif else dlam_all = replicate(3.9/2.35, ngoodlam)
  endelse
  fwhm_gal = 2.35*median(dlam_all)/(1.+science.zspec)

  ;;rebin lambda to have the same scale for all pixel
  ;mostly follow the example from ppxf_kinematics_example_sauron.pro
  lambda_in = science.lambda(goodlam)
  lambdiff = median(abs(ts_diff(lambda_in,1)))
  nspec    = fix((max(lambda_in)-min(lambda_in))/lambdiff)
  lambda   = findgen(nspec)*lambdiff+min(lambda_in)
  contdiv  = interpol(science.contdiv,science.lambda,lambda)
  contdivivar = science.contdivivar
  badpix = where(finite(contdivivar) eq 0 or contdivivar eq 0,cbadpix)
  if cbadpix gt 0 then contdivivar(badpix) = 1.e10
  noise = interpol(1./sqrt(science.contdivivar),science.lambda,lambda)
  ;switch to restframe
  lambda = lambda/(1.+science.zspec)
  lambrange = minmax(lambda)
  log_rebin,lambrange,contdiv,galaxy,loglam1,velscale=velscale
  noise  = noise/median(galaxy)
  galaxy = galaxy/median(galaxy)
  
  ;read ppxf's ssp library by Vazdekis1999
  vazdekis = file_search('/scr2/nichal/workspace2/idl/ppxf/spectra/Rbi1.30z*.fits',COUNT=nfiles)
  FWHM_tem = 1.8             ; Vazdekis spectra have a resolution FWHM of 1.8A.
  fits_read, vazdekis[0], ssp, h2
  lamRange2 = sxpar(h2,'CRVAL1') + [0d,sxpar(h2,'CDELT1')*(sxpar(h2,'NAXIS1')-1d)]
  log_rebin, lamRange2, ssp, sspNew, logLam2, VELSCALE=velScale
  templates = dblarr(n_elements(sspNew),nfiles)
  
  ;convolve the template with quadratic difference between galaxy and template
  FWHM_dif = SQRT(FWHM_gal^2 - FWHM_tem^2)
  sigma = FWHM_dif/2.355/sxpar(h2,'CDELT1') ; Sigma difference in pixels 

; IMPORTANT: To avoid spurious velocity offsets of the templates, the
; NPIXEL keyword in PSF_GAUSSIAN must be an odd integer as done below
;
  lsf = psf_Gaussian(NPIXEL=2*ceil(4*sigma)+1, ST_DEV=sigma, /NORM, NDIM=1)
  for j=0,nfiles-1 do begin
     fits_read, vazdekis[j], ssp
     ssp = convol(ssp,lsf)      ; Degrade template to SAURON resolution
    ; From IDL 8.1 one can use the following commented line instead of the 
    ; above one, and the line with PSF_GAUSSIAN is not needed any more.  
    ; ssp = gauss_smooth(ssp,sigma)    
     log_rebin, lamRange2, ssp, sspNew, VELSCALE=velScale
     templates[*,j] = sspNew/median(sspNew) ; Normalizes templates 
  endfor

  c = 299792.458d
  dv = (logLam2[0]-logLam1[0])*c ; km/s

  vel = 0d ; Initial estimate of the galaxy velocity in km/s
  goodPixels = ppxf_determine_goodPixels(logLam1,lamRange2,vel)

  start = [vel, 200d]           ; (km/s), starting guess for [V,sigma]
  ppxf, templates, galaxy, noise, velScale, start, sol, $
        GOODPIXELS=goodPixels, /PLOT, MOMENTS=4, DEGREE=4, $
        VSYST=dv, ERROR=error
  
  print, 'Formal errors:    dV    dsigma       dh3       dh4'
  print, error[0:3]*sqrt(sol[6]), FORMAT='(10x,2f10.1,2f10.3)'
  print, 'Best-fitting redshift z:', (science.zspec + 1)*(1 + sol[0]/c) - 1

  science.vdisp_ppxf = sol[1]
  science.vdisperr_ppxf = error[1]*sqrt(sol[6])
end

pro sps_fit::fit, science, noredraw=noredraw, nostatusbar=nostatusbar
    common sps_spec, sps, spsz, spsage
    common sps_iterproc, contiter
    common get_sps, dlam, dataivar, datalam, wonfit, contmask, normalize
    common toprint, agediff, zdiff
    common get_vdisp, dlamfit

    if ~keyword_set(nostatusbar) then widget_control, widget_info(self.base, find_by_uname='status'), set_value='Fitting ...'
    restlambda = science.lambda / (1d + science.zspec)
    znow = science.zspec
    reallambda = science.lambda
    nlambda = n_elements(reallambda)

    if min(science.dlam) gt 0.2 and max(science.dlam) lt 10.0 then begin
        dlam_all = science.dlam
    endif else begin
        specresfile = self.directory+'specres_poly.sav'
        if file_test(specresfile) then begin
            restore, specresfile
            dlam_all = poly(science.lambda/1000 - 7.8, specres_poly) / 2.35
        endif else dlam_all = replicate(3.9/2.35, nlambda)
    endelse

     ;FIT VELOCITY DISPERSION FIRST ;;;;;
    ;Use Mgb at 5150-5200 angstrom
    won = where(science.lambda/(1.+znow) ge 5150. and science.lambda/(1.+znow) lt 5200. and finite(science.contdiv) and finite(science.contdivivar), con)
    widget_control, widget_info(self.base, find_by_uname='2d'), get_value=index
    wset, index
    a = [znow,200.,0.2,0.2,0.2,1.] 
    xfit = science.lambda[won]
    yfit = science.contdiv[won]
    wfit = science.contdivivar[won] 
    dlamfit = dlam_all[won]
    fit = curvefit(xfit,yfit,wfit,a,sigmaa,chisq=chi2,function_name='get_vdisp',/noderivative)
    plot,xfit,yfit
    oplot,xfit,fit,color=fsc_color('red')
    science.vdisp_smm = a[1]
    science.vdisperr_smm = sigmaa[1]
    print, 'vdisp = ', science.vdisp_smm, science.vdisperr_smm
    print, 'redshift = ', a[0], sigmaa[0]
    
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ;FIT Whole Spectrum;;;;;;;;;;;;;
                                ;bestvdisp = wmean([science.vdisp_smm,science.vdisp_ppxf],[science.vdisperr_smm,science.vdisperr_ppxf])
                                ;vdisperror = sqrt(total([science.vdisperr_smm,science.vdisperr_ppxf]^2))
    bestvdisp = science.vdisp_smm
    vdisperror = science.vdisperr_smm
    pi = replicate({value:0d, fixed:0, limited:[1,1], limits:[0.D,0.D], parname:'', mpprint:0, mpformat:'', step:0d, tied:''}, 4)
    pi.value = double([-0.05, 3.0, bestvdisp,znow])
    pi[0].limits = minmax(spsz)
    pi[1].limits = minmax(spsage)
    pi[2].limits = [bestvdisp-2.*vdisperror,bestvdisp+2.*vdisperror]
    pi[3].limits = [-0.05,0.05]+znow
    pi.step = double([0.1, 0.5, 25.0,0.0002])
    pi.parname = ['    Z', '  age', 'vdisp','redshift']
    pi.mpformat = ['(D6.3)', '(D5.2)', '(D6.1)','(D6.3)']
   ; pi[1].limits[0] = 1.0


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
    nloop=0
    normalize=1
    bestchisq = 9999.
    bestvalue = [99.,99.,99.,99.]
    besterror = [99.,99.,99.,99.]
    maxnloop = 150
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, strtrim(science.objname, 2)+'  ('+strtrim(string(self.i+1, format='(I3)'), 2)+' / '+strtrim(string(self.nspec, format='(I3)'), 2)+')'
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, '  i Z/Z_sun   age sigma_v  redhift    chi^2  DOF   ZDIFF  AGEDIFF'
    print, '--- ------- ----- ------- ---------  -------- ---- -----  ------'
    curchisq = 0.
    while abs(zdiff) gt 0.0001 or abs(agediff) gt 0.0001 or abs(vdispdiff) gt 0.0001 or abs(redshfdiff) gt 0.0001 and nloop le maxnloop and curchisq lt 1.02*bestchisq do begin
        contiter++
        dlam = dlam_all
        dataivar = science.telldivivar*(median(science.telldiv))^2
        datalam = science.lambda
        wonfit = wontofit
        contmask = science.contmask
        pars = mpfitfun('get_sps_obs_082016', [2.*xmp[0]-xmp[1],xmp],[bestvdisp,ymp], [vdisperror,dymp], parinfo=pi, /nocatch, bestnorm=bestnorm, dof=dof, perror=perror, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=500, status=status, yfit=ympfit, iterproc='sps_iterproc')

        normalize=0
        zdiff = (pi[0].value-pars[0])/pi[0].value
        agediff = (pi[1].value-pars[1])/pi[1].value
        vdispdiff = (pi[2].value-pars[2])/pi[2].value
        redshfdiff = (pi[3].value-pars[3])/pi[2].value
        ;print,'z,age,vdisp diff:',zdiff,agediff,vdispdiff,nloop
        pi.value = pars

        nnan = where(finite(science.contdiv) and restlambda gt 3500 and restlambda lt 7400)
        wonfit  = nnan
        restlambda = science.lambda / (1d + pi[3].value)
        datalam = restlambda
        spsbestfitarr = get_sps_rest(restlambda, pi.value)                                
        spsbestfit=spsbestfitarr[*,1]

        
        spectoplot = science.telldiv/median(science.telldiv)
        !p.multi=[0,1,2]
        plot,restlambda,spsbestfitarr[*,1],/nodata,yrange=[0.1,max(spsbestfitarr[*,1])*1.2],xrange=[3500,6000]
        oplot,restlambda,spectoplot
        oplot,restlambda,spsbestfitarr[*,1],color=fsc_color('red')
       
        bkpt = slatec_splinefit(restlambda[won], science.contdiv[won] /spsbestfit[won], coeff, invvar=science.contdivivar[won], bkspace=150, upper=3, lower=3, /silent)
        if bkpt[0] eq -1 then begin
            pi.value = [-999d, -999d, -999d, -999d]
            perror = [-999d, -999d, -999d, -999d]
            science.spsspec = -999d
            science.spscont = -999d 
            break
        endif
        cont = slatec_bvalu(restlambda, bkpt, coeff)
        ympold = ymp
        ymp = science.contdiv[won] / cont[won]
        dymp = (science.contdivivar[won])^(-0.5) / cont[won]
        science.spscont = cont 
        nloop +=1
        plot,xmp/(1.+pi[3].value),ympold
        oplot,xmp/(1.+pi[3].value),ymp,color=fsc_Color('green')
        oplot,restlambda[won],spsbestfit[won],color=fsc_color('red')
        !p.multi=[0,1,1]

        if nloop eq maxnloop then print,'WARNING: MAX NLOOP REACHED!'

        curchisq = bestnorm/dof
        if curchisq lt bestchisq then begin
           bestchisq = curchisq
           bestvalue = pars
           besterror = perror
           bestspsbestfitarr = spsbestfitarr
           bestcont = cont
        endif
        
     endwhile
    print,agediff,zdiff,format='("--- ------- ----- ------- ---------  -------- ----",D9.6,2X,D9.6)'                            
    ;check if the last chisq is the best chisq
    if (curchisq-bestchisq)/bestchisq gt 0.005 then begin
       print,'THE WHILE LOOP HAS WALKED AWAY FROM THE BEST VALUES. BETTER CHECK YOUR PLOT'
       print,'The values used are:'
       print, bestvalue[0], bestvalue[1], bestvalue[2],bestvalue[3],bestchisq,format='(6X,D6.3,1X,D5.2,2X,D6.1,2x,D6.3,1X,D8.3)'
       science.goodfit = 1.
       spsbestfitarr = bestspsbestfitarr
       pi.value = bestvalue
       perror   = besterror
       science.spscont = bestcont
    endif

    science.nloop = nloop
    science.spsspec = spsbestfitarr[*,0] 
    science.spsspecfull = spsbestfitarr[*,1]
    science.spscontfull = spsbestfitarr[*,2]
    print, ' '
    
    done:
    science.feh = pi[0].value
    science.feherr = perror[0]
    science.age = pi[1].value
    science.ageerr = perror[1]
    science.zfit = pi[3].value
    science.zspec = pi[3].value

    science.vdisp = pi[2].value
    science.vdisperr = perror[2]
    ;calculate chisq
    science.chisq = total((spsbestfit[won]-science.contdiv[won]/science.spscont[won])^2*science.contdivivar[won]*(science.spscont[won])^2)/float(n_elements(won))

    self->statusbox, science=science
    if ~keyword_set(noredraw) then begin
        self->redraw
    endif
end

pro sps_fit::fitmg, science, noredraw=noredraw, nostatusbar=nostatusbar
    common sps_spec, sps, spsz, spsage
    common sps_iterproc, contiter
    common get_sps, dlam, dataivar, datalam, wonfit, contmask, normalize
    common toprint, agediff, zdiff

    if ~keyword_set(nostatusbar) then widget_control, widget_info(self.base, find_by_uname='status'), set_value='Fitting Mg...'
    restlambda = science.lambda / (1d + science.zspec)
    znow = science.zfit
    reallambda = science.lambda
    nlambda = n_elements(reallambda)
    
    if min(science.dlam) gt 0.2 and max(science.dlam,/nan) lt 10.0 then begin
        dlam_all = science.dlam
    endif else begin
        specresfile = self.directory+'specres_poly.sav'
        if file_test(specresfile) then begin
            restore, specresfile
            dlam_all = poly(science.lambda/1000 - 7.8, specres_poly) / 2.35
        endif else dlam_all = replicate(3.9/2.35, nlambda)
    endelse
    
    pi = replicate({value:0d, fixed:0, limited:[1,1], limits:[0.D,0.D], parname:'', mpprint:0, mpformat:'', step:0d, tied:''}, 4)
    pi.value = double([science.feh, science.age, science.vdisp,science.zfit])
    pi[0].limits = minmax(spsz)
    pi[1].limits = minmax(spsage)
    pi[2].limits = [0.0, 1000.0]
    pi[3].limits = [-0.05,0.05]+znow
    pi.step = double([0.1, 0.5, 25.0,0.0002])
    pi.parname = ['    Z', '  age', 'vdisp','redshift']
    pi.mpformat = ['(D6.3)', '(D5.2)', '(D6.1)','(D6.3)']
    pi.fixed = [0,1,1,1]

    won = where(science.fitmgmask eq 1 and finite(science.contdiv) and finite(science.contdivivar) and science.contdivivar gt 0, con)
    if con lt 10 then begin
        pi.value = [-999d, -999d, -999d, -999d]
        perror = [-999d, -999d, -999d, -999d]
        science.spsspec = -999d
        goto, done
    endif
    
    xmp = reallambda[won]
    ymp = science.contdiv[won]/science.spscont[won]
    dymp = (science.contdivivar[won])^(-0.5)/science.spscont[won]
    wontofit = won
    contiter = 0
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, strtrim(science.objname, 2)+'  ('+strtrim(string(self.i+1, format='(I3)'), 2)+' / '+strtrim(string(self.nspec, format='(I3)'), 2)+') Mg'
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, '  i Z/Z_sun   age sigma_v  redhift    chi^2  DOF'
    print, '--- ------- ----- ------- ---------  -------- ----'
    nloop = 0
    normalize = 0

    contiter++
    dlam = dlam_all
    dataivar = science.telldivivar*(median(science.telldiv))^2
    datalam  = science.lambda
    wonfit   = wontofit
    contmask = science.contmask
    pars = mpfitfun('get_sps_obs', xmp, ymp, dymp, parinfo=pi, /nocatch, bestnorm=bestnorm, dof=dof, perror=perror, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=1000, status=status, yfit=ympfit, iterproc='sps_iterproc')
        
    pi.value = pars
    restlambda = science.lambda / (1d + pi[3].value)
    datalam = restlambda
    spsbestfitarr = get_sps_rest(restlambda, pi.value)
 
    science.spsspecmg = spsbestfitarr[*,0]
    science.spsspecfullmg = spsbestfitarr[*,1]
    science.spscontfullmg = spsbestfitarr[*,2]
    print, ' '
    
    science.mg = pi[0].value
    science.mgerr = perror[0]
 
    ;calculate chisq
    science.chisqmg = total((science.spsspecfullmg[won]-science.contdiv[won]/science.spscont[won])^2*science.contdivivar[won]*(science.spscont[won])^2)/float(n_elements(won))
    done:
    self->statusbox, science=science
    if ~keyword_set(noredraw) then begin
        self->redraw
     endif

 end

pro sps_fit::fitfe, science, noredraw=noredraw, nostatusbar=nostatusbar
    common sps_spec, sps, spsz, spsage
    common sps_iterproc, contiter
    common get_sps, dlam, dataivar, datalam, wonfit, contmask, normalize
    common toprint, agediff, zdiff
    if ~keyword_set(nostatusbar) then widget_control, widget_info(self.base, find_by_uname='status'), set_value='Fitting Fe...'
    restlambda = science.lambda / (1d + science.zspec)
    znow = science.zfit
    reallambda = science.lambda
    nlambda = n_elements(reallambda)
    
    if min(science.dlam) gt 0.2 and max(science.dlam,/nan) lt 10.0 then begin
        dlam_all = science.dlam
    endif else begin
        specresfile = self.directory+'specres_poly.sav'
        if file_test(specresfile) then begin
            restore, specresfile
            dlam_all = poly(science.lambda/1000 - 7.8, specres_poly) / 2.35
        endif else dlam_all = replicate(3.9/2.35, nlambda)
    endelse
    
    pi = replicate({value:0d, fixed:0, limited:[1,1], limits:[0.D,0.D], parname:'', mpprint:0, mpformat:'', step:0d, tied:''}, 4)
    pi.value = double([science.feh, science.age, science.vdisp,science.zfit])
    pi[0].limits = minmax(spsz)
    pi[1].limits = minmax(spsage)
    pi[2].limits = [0.0, 1000.0]
    pi[3].limits = [-0.05,0.05]+znow
    pi.step = double([0.1, 0.5, 25.0,0.0002])
    pi.parname = ['    Z', '  age', 'vdisp','redshift']
    pi.mpformat = ['(D6.3)', '(D5.2)', '(D6.1)','(D6.3)']
    pi.fixed = [0,1,1,1]

    won = where(science.fitfemask eq 1 and finite(science.contdiv) and finite(science.contdivivar) and science.contdivivar gt 0, con)
    if con lt 10 then begin
        pi.value = [-999d, -999d, -999d, -999d]
        perror = [-999d, -999d, -999d, -999d]
        science.spsspec = -999d
        goto, done
    endif
    
    xmp = reallambda[won]
    ymp = science.contdiv[won]/science.spscont[won]
    dymp = (science.contdivivar[won])^(-0.5)/science.spscont[won]
    wontofit = won
    zdiff = 1.0
    agediff = 1.0
    vdispdiff = 1.0
    redshfdiff = 1.0
    contiter = 0
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, strtrim(science.objname, 2)+'  ('+strtrim(string(self.i+1, format='(I3)'), 2)+' / '+strtrim(string(self.nspec, format='(I3)'), 2)+') Fe'
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, '  i Z/Z_sun   age sigma_v  redhift    chi^2  DOF'
    print, '--- ------- ----- ------- ---------  -------- ----'
    nloop = 0
    normalize = 0

    contiter++
    dlam = dlam_all
    dataivar = science.telldivivar*(median(science.telldiv))^2
    datalam  = science.lambda
    wonfit   = wontofit
    contmask = science.contmask
    pars = mpfitfun('get_sps_obs',xmp, ymp, dymp, parinfo=pi, /nocatch, bestnorm=bestnorm, dof=dof, perror=perror, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=1000, status=status, yfit=ympfit, iterproc='sps_iterproc')
        
    pi.value = pars
    restlambda = science.lambda / (1d + pi[3].value)
    datalam = restlambda
    
    spsbestfitarr = get_sps_rest(restlambda, pi.value)
 
    science.spsspecfe = spsbestfitarr[*,0]
    science.spsspecfullfe = spsbestfitarr[*,1]
    science.spscontfullfe = spsbestfitarr[*,2]
    print, ' '
    
    science.fe = pi[0].value
    science.feerr = perror[0]
 
    ;calculate chisq
    science.chisqfe = total((science.spsspecfullfe[won]-science.contdiv[won]/science.spscont[won])^2*science.contdivivar[won]*(science.spscont[won])^2)/float(n_elements(won))

    done:
    self->statusbox, science=science
    if ~keyword_set(noredraw) then begin
        self->redraw
    endif
 end

pro sps_fit::fit_all
    scienceall = *self.science
    curi = self.i
    for i=0,self.nspec-1 do begin
        self.i = i
        self->default_range
        science = scienceall[self.i]
        if science.good eq 0 then continue
        self->vdisp_ppxf,science
        scienceall[self.i] = science
        self->fit, science
        scienceall[self.i] = science
        self->maskfe, science
        self->fitfe, science
        scienceall[self.i] = science
        self->maskmg, science
        self->fitmg, science
        scienceall[self.i] = science
    endfor
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    self->writescience
    self.i = curi
    science = scienceall[self.i]
    self->statusbox, science=science
    self->redraw    
end


; =================
pro sps_fit_event, ev
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_uvalue=value
    if (n_elements(value) eq 0) then value = ''
    name = strmid(tag_names(ev, /structure_name), 7, 4)

    case (name) of
        'BUTT': obj->handle_button, ev
        'TEXT': obj->handle_text, ev
        'TRAC': obj->handle_tracking, ev
        'COMB': obj->handle_combobox, ev
        'DRAW': begin
            if ev.type eq 0 then obj->handle_draw_click, ev
            if ev.type eq 2 then obj->handle_motion, ev
            if ev.type eq 5 then obj->handle_draw_key, ev
        end
        'DONE': widget_control, ev.top, /destroy
        else: begin
            case (value) of
                'good': obj->toggle_good
                else: obj->redraw
            endcase
            widget_control, widget_info(ev.top, find_by_uname='spec'), /input_focus
        end
    endcase
end


; ================ BUTTONS ================
function line_event, ev
    self->redraw
end


pro sps_fit::handle_button, ev
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_uvalue=uvalue
    
    case (uvalue) of
        'back': self->newspec, increment=-1
        'next': self->newspec, increment=1
        'default_range': begin
            self->default_range
            self->redraw
        end
        'reprepare': self->reprepare
        'indices': begin
            scienceall = *self.science
            science = scienceall[self.i]
            self->indices, science
            scienceall[self.i] = science
            ptr_free, self.science
            self.science = ptr_new(scienceall)
            self->statusbox
            ;self->redraw
         end
        'fit': begin
           scienceall = *self.science
           science = scienceall[self.i]
           self->fit, science, /noredraw
           scienceall[self.i] = science
           ptr_free, self.science
           self.science = ptr_new(scienceall)
           self->redraw
        end
        'fit[Mg/Fe]': begin
           scienceall = *self.science
           science = scienceall[self.i]
           self->maskfe,science
           self->fitfe, science, /noredraw
           scienceall[self.i] = science
           self->maskmg,science
           self->fitmg, science, /noredraw
           scienceall[self.i] = science
           ptr_free, self.science
           self.science = ptr_new(scienceall)
           self->redraw
        end
        'reprepare_all': self->reprepare_all
        'fit_all': self->fit_all
        'default_cont': self->default_cont
        'default_mask': self->default_mask
        'default_maskall': self->default_maskall
        'default_goodspec':self->default_goodspec
        'backward': self->step, -1
        'forward': self->step, 1
        'good': self->toggle_good
        'blue': self->lambdarange, /blue
        'red': self->lambdarange, /red
        'save': self->writescience
        'exit': begin
            self->writescience
            widget_control, ev.top, /destroy
        end
        else:
    endcase
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end


pro sps_fit::step, increment
    nmodes = 5
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode eq 1 and increment gt 0: mode = 3
        mode eq 3 and increment gt 0: mode = 4
        mode eq 4 and increment lt 0: mode = 3
        mode eq 3 and increment lt 0: mode = 1
        else: mode += increment
    endcase
    curi = self.i
    if mode eq -1 then begin
        self->newspec, increment=-1, /noredraw
        if self.i eq curi then return
        mode = nmodes-1
        ;self->default_range
    endif
    if mode eq nmodes then begin
        self->newspec, increment=1, /noredraw
        if self.i eq curi then return
        mode = 1
        ;self->default_range
    endif
    widget_control, widget_info(self.base, find_by_uname='mode'), set_value=mode
    self.keystate = 0
    self->redraw
end


pro sps_fit::toggle_good
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Updating database ...'
    scienceall = *self.science
    widget_control, widget_info(self.base, find_by_uname='good'), get_value=good
    scienceall[self.i].goodsky = good[0]
    scienceall[self.i].good = good[1]
    scienceall[self.i].goodfit = good[2]
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end


pro sps_fit::output_objname
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Outputting object name to output.dat ...'
    comment = textbox(title='Provide a Comment (optional)', group_leader=self.base, label='Comment: ', cancel=cancelled, xsize=59, value='')
    if ~cancelled then begin
        science = (*self.science)[self.i]
        openw, lun, 'output.dat', /append, /get_lun
        printf, lun, science.objname, comment, format='(A-20,1X,A-59)'
        close, lun
        free_lun, lun
    endif
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end


pro sps_fit::newspec, increment=increment, noredraw=noredraw
    newi = self.i
    ;good = 0
    ;while(~good) do begin
        newi += increment
        if newi lt 0 then begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='This is the first spectrum.'
            good = 1
            return
        endif
        if newi gt self.nspec-1 then begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='This is the last spectrum.'
            good = 1
            return
        endif
    ;    good = (*self.science)[newi].good eq 1
    ;endwhile
    widget_control, widget_info(self.base, find_by_uname='filelist'), set_combobox_select=newi
    self.i = newi
    self->statusbox
    self.ylim = minmax((*self.science)[self.i].spec)
    self->skyylim
    ;self->default_range
    
    self.keystate = 0
    if ~keyword_set(noredraw) then self->redraw

    common slit, bgood, rgood, bslit, rslit, bwidth, rwidth, bheight, rheight
    bgood = 0
    rgood = 0
    science = (*self.science)[self.i]
    if strmid(science.mask,0,4) eq '0024' then maskname = '0024' else maskname = science.mask
    if strmid(science.mask,0,5) eq '0024b' then maskname = science.mask
    bslitfile = file_dirname(science.spec1dfile, /mark_directory)+'slit.'+maskname+'.'+string(science.slit, format='(I03)')+'B.fits.gz'
    rslitfile = file_dirname(science.spec1dfile, /mark_directory)+'slit.'+maskname+'.'+string(science.slit, format='(I03)')+'R.fits.gz'
    widget_control, widget_info(self.base, find_by_uname='2d'), get_value=index
    wset, index
    loadct, 0, /silent
    tv, dblarr(1600, 300), 0, 0
    if file_test(bslitfile) then bslit = mrdfits(bslitfile, 1, /silent)
    if file_test(rslitfile) then rslit = mrdfits(rslitfile, 1, /silent)
    case 1 of
        file_test(bslitfile) and file_test(rslitfile): begin
            bflux = bslit.flux
            rflux = rslit.flux
            width = min([(size(bflux, /dimensions))[1], (size(rflux, /dimensions))[1]])
            flux = [bflux[*,0:width-1], rflux[*,0:width-1]]
        end
        file_test(bslitfile): flux = bslit.flux
        file_test(rslitfile): flux = rslit.flux
        else: return
    endcase
    med = median(flux)
    sig = stddev(flux)
    max_value = (med + (10 * sig)) < max(flux)
    min_value = (med - (2 * sig))  > min(flux)
    if (finite(min_value) EQ 0) then min_value = image_min
    if (finite(max_value) EQ 0) then max_value = image_max
    if (min_value GE max_value) then begin
        min_value = min_value - 1
        max_value = max_value + 1
    endif

    if file_test(bslitfile) then begin
        bflux = bytscl(bslit.flux, /nan, min=min_value, max=max_value, top=255) + 8
        bflux = -1.0*bflux + 255.0
        bheight = 1600./2047.*(size(bslit.dlambda, /dimensions))[1]
        if bheight gt 74 then begin
            bwidth = round(74./bheight*1600.)
            bheight = 74
        endif else begin
            bwidth = 1600
            bheight = round(bheight)
        endelse
        bflux1 = congrid(bflux[0:2047,*], bwidth, bheight, /interp)
        bflux2 = congrid(bflux[2048:4095,*], bwidth, bheight, /interp)
        tv, bflux1, 0, 225
        tv, bflux2, 0, 150
        bgood = 1
    endif else bgood = 0
    if file_test(rslitfile) then begin
        rflux = bytscl(rslit.flux, /nan, min=min_value, max=max_value, top=247) + 8
        rflux = -1.0*rflux + 255.0
        rheight = 1600./2047.*(size(rslit.dlambda, /dimensions))[1]
        if rheight gt 74 then begin
            rwidth = round(74./rheight*1600.)
            rheight = 74
        endif else begin
            rwidth = 1600
            rheight = round(rheight)
        endelse
        rflux1 = congrid(rflux[0:2047,*], rwidth, rheight, /interp)
        rflux2 = congrid(rflux[2048:4095,*], rwidth, rheight, /interp)
        tv, rflux1, 0, 75
        tv, rflux2, 0, 0
        rgood = 1
    endif else rgood = 0

    if bgood eq 1 and rgood eq 1 then begin
        n = n_elements(*self.tellstart)
        for i=0,n-1 do begin
            x1 = interpol(dindgen(4096), bslit.lambda0, (*self.tellstart)[i])
            offset1 = 0
            case 1 of
                x1 lt 0: x1 = -1
                x1 le 2047: offset1 = 225
                x1 le 4095: begin
                    offset1 = 150
                    x1 -= 2047
                end
                x1 gt 4095: begin
                    x1 = interpol(dindgen(4096), rslit.lambda0, (*self.tellstart)[i])
                    case 1 of
                        x1 lt 0: x1 = -1
                        x1 le 2047: offset1 = 75
                        x1 le 4095: begin
                            offset1 = 0
                            x1 -= 2047
                        end
                        x1 gt 4095: x1 = -1
                    endcase
                end
            endcase
            if x1 gt 0 then begin
                width1 = offset1 ge 150 ? bwidth : rwidth
                x1 *= width1 / 2047.
                height1 = offset1 ge 150 ? bheight : rheight
            endif else begin
                width1 = 0
                height1 = 0
            endelse

            x2 = interpol(dindgen(4096), bslit.lambda0, (*self.tellend)[i])
            offset2 = 0
            case 1 of
                x2 lt 0: begin
                    x2 = 2047
                    offset2 = offset1
                end
                x2 le 2047: offset2 = 225
                x2 le 4095: begin
                    offset2 = 150
                    x2 -= 2047
                end
                x2 gt 4095: begin
                    x2 = interpol(dindgen(4096), rslit.lambda0, (*self.tellend)[i])
                    case 1 of
                        x2 lt 0: begin
                            x2 = 2047
                            offset2 = offset1
                        end
                        x2 le 2047: offset2 = 75
                        x2 le 4095: begin
                            offset2 = 0
                            x2 -= 2047
                        end
                        x2 gt 4095: x2 = 2047
                    endcase
                end
            endcase
            if x2 gt 0 then begin
                width2 = offset2 ge 150 ? bwidth : rwidth
                x2 *= width2 / 2047.
                height2 = offset2 ge 150 ? bheight : rheight
            endif else begin
                width2 = 0
                height2 = 0
            endelse

            case 1 of
                x1 lt 0 and x2 lt 0: break
                x1 ge width1 and x2 ge width2: break
                offset1 eq offset2: begin
                    plots, [x1 > 0, x2 < width2], [0, 0]+height1-3+offset1, color=fsc_color('green'), /device, thick=(*self.tellthick)[i]
                    plots, [x1 > 0, x2 < width2], [0, 0]+2+offset1, color=fsc_color('green'), /device, thick=(*self.tellthick)[i]
                end
                offset1 gt offset2: begin
                    plots, [x1 > 0, width2], [0, 0]+height1-3+offset1, color=fsc_color('green'), /device, thick=(*self.tellthick)[i]
                    plots, [x1 > 0, width2], [0, 0]+2+offset1, color=fsc_color('green'), /device, thick=(*self.tellthick)[i]
                    plots, [0, x2 < width2], [0, 0]+height2-3+offset2, color=fsc_color('green'), /device, thick=(*self.tellthick)[i]
                    plots, [0, x2 < width2], [0, 0]+2+offset2, color=fsc_color('green'), /device, thick=(*self.tellthick)[i]
                end
                else:
            endcase
        endfor

        n = n_elements(*self.linewaves)
        for i=0,n-1 do begin
            x = interpol(dindgen(4096), bslit.lambda0, (*self.linewaves)[i] * (1d + science.zspec))
            case 1 of
                x lt 0: x = -1
                x le 2047: offset = 225
                x le 4095: begin
                    offset = 150
                    x -= 2047
                end
                x gt 4095: begin
                    x = interpol(dindgen(4096), rslit.lambda0, (*self.linewaves)[i] * (1d + science.zspec))
                    case 1 of
                        x lt 0: x = -1
                        x le 2047: offset = 75
                        x le 4095: begin
                            offset = 0
                            x -= 2047
                        end
                        x gt 4095: x = -1
                    endcase
                end
            endcase
            if x lt 0 then continue
            x *= (offset ge 150 ? bwidth : rwidth) / 2047.
            height = offset ge 150 ? bheight : rheight
            plots, [x, x], [0, 5]+offset, color=fsc_color((*self.linecolors)[i]), /device
            plots, [x, x], [height-6, height-1]+offset, color=fsc_color((*self.linecolors)[i]), /device
            xyouts, x+3, height-7+offset, (*self.linenames)[i], color=fsc_color((*self.linecolors)[i]), orientation=90, alignment=1, /device
        endfor
     endif
end


pro sps_fit::skyylim
    wsky = where((*self.science)[self.i].skylinemask ne -1, csky)
    if csky gt 0 then begin
        skyfit = (*self.science)[self.i].skyfit[wsky,*]
        yminmax = (*self.science)[self.i].dlam
        if skyfit[0,0] ne -1 then yminmax = [yminmax, skyfit[wsky,1]+skyfit[wsky,2], skyfit[wsky,1]-skyfit[wsky,2]]    
        self.skyylim = minmax(yminmax*2.35)    
    endif else begin
        self.skyylim = [0.9, 1.7]
    endelse
end


pro sps_fit::default_range, update=update
    if ~keyword_set(update) then begin
        self.ylim = minmax((*self.science)[self.i].spec)
        self.ylim[1] *= 1.1
        self.divylim = [-1.0, 2.5]
        self->skyylim
        self.lambdalim = (minmax((*self.science)[self.i].lambda / (1d + (*self.science)[self.i].zspec)) < 9100) > 2000
        self.lambdalim[0] >= 2000.
        self.lambdalim[1] <= 8938. / (1d + (*self.science)[self.i].z)
    endif
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode le 1: begin
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.ylim[0], format='(g8.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.ylim[1], format='(g8.2)'), /rem)
        end
        mode eq 3: begin
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.skyylim[0], format='(D5.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.skyylim[1], format='(D5.2)'), /rem)
        end
        else: begin
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.divylim[0], format='(D5.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.divylim[1], format='(D5.2)'), /rem)
        end
    endcase
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    zl = mode lt 4 ? (*self.science)[self.i].z : 0.0
    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
end


pro sps_fit::lambdarange, red=red, blue=blue
    if keyword_set(red)+keyword_set(blue) ne 1 then message, 'You must specify red or blue.'
    lrange = self.lambdalim[1] - self.lambdalim[0]
    if keyword_set(blue) then begin
        if self.lambdalim[0] lt min((*self.science)[self.i].lambda) then begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='The is bluest part of the spectrum.'
            return
        endif
        llownew = self.lambdalim[0] - 0.6*lrange
        lhighnew = self.lambdalim[1] - 0.6*lrange
    endif
    if keyword_set(red) then begin
        if self.lambdalim[1] gt max((*self.science)[self.i].lambda) then begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='The is reddest part of the spectrum.'
            return
        endif
        llownew = self.lambdalim[0] + 0.6*lrange
        lhighnew = self.lambdalim[1] + 0.6*lrange
    endif
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    zl = mode lt 5 ? (*self.science)[self.i].z : 0.0
    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
    self->lambdalow, llownew, /noredraw
    self->lambdahigh, lhighnew
end


; ============== TEXT BOXES =============
pro sps_fit::handle_text, ev
    widget_control, ev.id, get_uvalue=uvalue
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_value=val

    case (uvalue) of
        'ylow': self->ylow, val
        'yhigh': self->yhigh, val
        'lambdalow': self->lambdalow, val
        'lambdahigh': self->lambdahigh, val
        else: 
    end
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end

pro sps_fit::lambdalow, lambdalow, noredraw=noredraw
    ;if lambdalow gt max((*self.science)[self.i].lambda) then begin
    ;    lambdalow = max((*self.science)[self.i].lambda) - (self.lambdalim[1] - lambdalow)
    ;    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    ;    zl = mode lt 4 ? (*self.science)[self.i].z : 0.0
    ;    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
    ;endif
    self.lambdalim[0] = lambdalow
    if ~keyword_set(noredraw) then self->redraw
end


pro sps_fit::lambdahigh, lambdahigh, noredraw=noredraw
    ;if lambdahigh lt min((*self.science)[self.i].lambda) then begin
    ;    lambdahigh = min((*self.science)[self.i].lambda) + (lambdahigh - self.lambdalim[0])
    ;    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    ;    zl = mode lt 4 ? (*self.science)[self.i].z : 0.0
    ;    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
    ;endif
    self.lambdalim[1] = lambdahigh
    if ~keyword_set(noredraw) then self->redraw
end


pro sps_fit::ylow, ylow, noredraw=noredraw
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode le 1: self.ylim[0] = ylow
        mode eq 3: self.skyylim[0] = ylow
        else: self.divylim[0] = ylow
    endcase
    if ~keyword_set(noredraw) then self->redraw
end


pro sps_fit::yhigh, yhigh
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode le 1: self.ylim[1] = yhigh
        mode eq 3: self.skyylim[1] = yhigh
        else: self.divylim[1] = yhigh
    endcase
    self->redraw
end


pro sps_fit::smooth, smooth    
    leftover = smooth mod 2
    smooth = fix(smooth)
    if leftover ge 1 then smooth += 1
    smooth >= 0
    self.smooth = smooth
    widget_control, widget_info(self.base, find_by_uname='smooth'), set_value=strcompress(self.smooth, /rem)
    self->redraw
end


; =============== TRACKING ==============
pro sps_fit::handle_tracking, ev
    uname = widget_info(ev.id, /uname)

    case (uname) of
        '2d': begin
            if ev.enter eq 0 then begin
                widget_control, widget_info(self.base, find_by_uname='obswave'), set_value='     '
                widget_control, widget_info(self.base, find_by_uname='restwave'), set_value='     '
                widget_control, widget_info(self.base, find_by_uname='yslit'), set_value='     '
            endif
        end
        else: 
    end
end


; ================ MOTION ===============
pro sps_fit::handle_motion, ev
    uname = widget_info(ev.id, /uname)

    case (uname) of
        '2d': begin
            common slit, bgood, rgood, bslit, rslit, bwidth, rwidth, bheight, rheight
            if bgood eq 0 or rgood eq 0 then return
            valid = 1
            case 1 of
                ev.y ge 225 and ev.y lt 225+bheight and ev.x lt bwidth: begin
                    slit = bslit
                    bheightorig = (size(bslit.dlambda, /dimensions))[1]
                    y = double(ev.y)-225.
                    if bheightorig gt 74 then y *= double(bheightorig) / 74.
                    x = double(ev.x) * 2047. / double(bwidth)
                end
                ev.y ge 150 and ev.y lt 150+bheight and ev.x lt bwidth: begin
                    slit = bslit
                    bheightorig = (size(bslit.dlambda, /dimensions))[1]
                    y = double(ev.y)-150.
                    if bheightorig gt 74 then y *= double(bheightorig) / 74.
                    x = (double(ev.x) * 2047. / double(bwidth)) + 2047.
                end
                ev.y ge 75 and ev.y lt 75+rheight and ev.x lt rwidth: begin
                    slit = rslit
                    rheightorig = (size(rslit.dlambda, /dimensions))[1]
                    y = double(ev.y)-75.
                    if rheightorig gt 74 then y *= double(rheightorig) / 74.
                    x = double(ev.x) * 2047. / double(rwidth)
                end
                ev.y ge 0 and ev.y lt rheight and ev.x lt rwidth: begin
                    slit = rslit
                    rheightorig = (size(rslit.dlambda, /dimensions))[1]
                    y = double(ev.y)
                    if rheightorig gt 74 then y *= double(rheightorig) / 74.
                    x = (double(ev.x) * 2047. / double(rwidth)) + 2047.
                end
                else: valid = 0
            endcase
            if ~valid then begin
                widget_control, widget_info(self.base, find_by_uname='obswave'), set_value='     '
                widget_control, widget_info(self.base, find_by_uname='restwave'), set_value='     '
                widget_control, widget_info(self.base, find_by_uname='yslit'), set_value='     '
            endif else begin
                lambda = slit.lambda0[x] + slit.dlambda[x, y]
                widget_control, widget_info(self.base, find_by_uname='obswave'), set_value='obs = '+string(lambda, format='(D7.1)')
                widget_control, widget_info(self.base, find_by_uname='restwave'), set_value='rest = '+string(lambda / (1.0 + (*self.science)[self.i].z), format='(D7.1)')
                widget_control, widget_info(self.base, find_by_uname='yslit'), set_value='y = '+string(y, format='(I3)')
            endelse
        end
        else: 
    end
end


; ============== COMBOBOX  =============
pro sps_fit::handle_combobox, ev
    ;widget_control, ev.id, get_uvalue=uvalue
    ;widget_control, ev.top, get_uvalue=obj
    self.i = ev.index
    self->statusbox
    self.ylim = minmax((*self.science)[self.i].spec)
    self->skyylim
    self.keystate = 0
    self->redraw
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end


; ============= DRAW CLICK =============
pro sps_fit::handle_draw_click, ev
    click_coords = convert_coord(ev.x, ev.y, /device, /to_data)
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    if mode lt 4 then click_coords /= 1d + (*self.science)[self.i].z
    case ev.modifiers of
        0: begin
            case ev.press of
                1: lrange = (self.lambdalim[1] - self.lambdalim[0]) / 2.
                4: lrange = (self.lambdalim[1] - self.lambdalim[0]) * 2.
                else: begin
                    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
                    return
                end
            endcase
            llownew = click_coords[0] - lrange/2.
            lhighnew = click_coords[0] + lrange/2.
            zl = mode lt 4 ? (*self.science)[self.i].z : 0.0
            widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
            self->lambdalow, llownew, /noredraw
            self->lambdahigh, lhighnew
        end
        1: begin
            if ev.press ne 1 then begin
                widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
                return
            endif
            lrange = self.lambdalim[1] - self.lambdalim[0]
            llownew = click_coords[0] - lrange/2.
            lhighnew = click_coords[0] + lrange/2.        
            zl = mode lt 4 ? (*self.science)[self.i].z : 0.0
            widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
            self->lambdalow, llownew, /noredraw
            self->lambdahigh, lhighnew
        end
        2: begin
            widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
            case mode of
                mode le 1: ylim = self.ylim
                mode eq 3: ylim = self.skyylim
                else: ylim = self.divylim
            endcase
            case ev.press of
                1: yrange = (ylim[1] - ylim[0]) / 2.
                4: yrange = (ylim[1] - ylim[0]) * 2.
                else: begin
                    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
                    return
                end
            endcase
            ylownew = click_coords[1] - yrange/2.
            yhighnew = click_coords[1] + yrange/2.
            case 1 of
                mode le 1: begin
                    widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(ylownew, format='(g8.2)'), /rem)
                    widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(yhighnew, format='(g8.1)'), /rem)
                end
                else: begin
                    widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(ylownew, format='(D5.2)'), /rem)
                    widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(yhighnew, format='(D5.1)'), /rem)
                end
            endcase
            self->ylow, ylownew, /noredraw
            self->yhigh, yhighnew
        end
        else: begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
            return            
        end
    endcase
end


; ============= DRAW KEYS ==============
pro sps_fit::handle_draw_key, ev
    if ev.press ne 1 then return
    key = string(ev.ch)
    coords = convert_coord(ev.x, ev.y, /device, /to_data)
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    science = (*self.science)[self.i]
    case mode of
        -1: begin        ;smoothed
            case key of
                'b': self->newspec, increment=-1
                'n': self->newspec, increment=1
                'g': begin
                    widget_control, widget_info(self.base, find_by_uname='good'), set_value=[science.goodsky, ~science.good]
                    self->toggle_good
                end
                'f': begin
                    self->default_range
                    self->redraw
                end
                'o': self->output_objname
                else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
            endcase
        end

        0: begin        ;telluric cross-correlation
            case key of
                'b': self->newspec, increment=-1
                'n': self->newspec, increment=1
                'g': begin
                    widget_control, widget_info(self.base, find_by_uname='good'), set_value=[science.goodsky, ~science.good]
                    self->toggle_good
                end
                'f': begin
                    self->default_range
                    self->redraw
                end
                'o': self->output_objname
                else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
            endcase
        end

        1: begin        ;continuum fit
            case key of
                'g': begin
                    case self.keystate of
                        0: begin
                            widget_control, widget_info(self.base, find_by_uname='good'), set_value=[science.goodsky, ~science.good]
                            self->toggle_good
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'b': begin
                    case self.keystate of
                        0: begin
                            self->newspec, increment=-1
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'n': begin
                    case self.keystate of
                        0: begin
                            self->newspec, increment=1
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'f': begin
                    case self.keystate of
                        0: begin
                            self->default_range
                            self->redraw
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'q': begin
                    case 1 of
                        self.keystate eq 1 or self.keystate eq 2: begin
                            self.keystate = 0
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Continuum mask modification cancelled.'
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'e': begin
                    case self.keystate of
                        0: begin
                            self.lambda1 = coords[0]
                            self.keystate = 1
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value="Press 'e' again to exclude wavelength region from continuum mask."
                        end
                        1: begin
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Updating database ...'
                            lambda1 = self.lambda1 < coords[0]
                            lambda2 = self.lambda1 > coords[0]
                            self.keystate = 0
                            scienceall = *self.science                            
                            w = where(scienceall[self.i].lambda gt lambda1 and scienceall[self.i].lambda lt lambda2, c)
                            if c eq 0 then begin
                                widget_control, widget_info(self.base, find_by_uname='status'), set_value='This wavelength region is invalid.'
                            endif else begin                            
                                scienceall[self.i].contmask[w] = 0
                                ptr_free, self.science
                                self.science = ptr_new(scienceall)
                                self->redraw
                            endelse
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'i': begin
                    case self.keystate of
                        0: begin
                            self.lambda1 = coords[0]
                            self.keystate = 2
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value="Press 'i' again to include wavelength region in continuum mask."
                        end
                        2: begin
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Updating database ...'
                            lambda1 = self.lambda1 < coords[0]
                            lambda2 = self.lambda1 > coords[0]
                            self.keystate = 0
                            scienceall = *self.science                            
                            w = where(scienceall[self.i].lambda gt lambda1 and scienceall[self.i].lambda lt lambda2, c)
                            if c eq 0 then begin
                                widget_control, widget_info(self.base, find_by_uname='status'), set_value='This wavelength region is invalid.'
                            endif else begin                            
                                scienceall[self.i].contmask[w] = 1
                                ptr_free, self.science
                                self.science = ptr_new(scienceall)
                                self->redraw
                            endelse
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase

                end
                'o': self->output_objname
                else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
            endcase
        end

        2: begin        ;continuum division
            case key of
                'g': begin
                    widget_control, widget_info(self.base, find_by_uname='good'), set_value=[science.goodsky, ~science.good]
                    self->toggle_good
                end
                'b': self->newspec, increment=-1
                'n': self->newspec, increment=1
                'f': begin
                    self->default_range
                    self->redraw
                end
                'o': self->output_objname
                else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
            endcase
        end

        3: begin        ;sky line fit
            case key of
                's': begin
                    widget_control, widget_info(self.base, find_by_uname='good'), set_value=[~science.goodsky, science.good]
                    self->toggle_good
                end
                'g': begin
                    widget_control, widget_info(self.base, find_by_uname='good'), set_value=[science.goodsky, ~science.good]
                    self->toggle_good
                end
                'b': self->newspec, increment=-1
                'n': self->newspec, increment=1
                'a': self->add_skyline, coords[0]
                'd': self->delete_skyline, coords[0]
                'f': begin
                    self->default_range
                    self->redraw
                end
                'o': self->output_objname
                else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
            endcase
        end

        4: begin        ;pixel mask
            case key of
                'g': begin
                    case self.keystate of
                        0: begin
                            widget_control, widget_info(self.base, find_by_uname='good'), set_value=[science.goodsky, ~science.good]
                            self->toggle_good
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'b': begin
                    case self.keystate of
                        0: begin
                            self->newspec, increment=-1
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'n': begin
                    case self.keystate of
                        0: begin
                            self->newspec, increment=1
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'f': begin
                    case self.keystate of
                        0: begin
                            self->default_range
                            self->redraw
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'h': begin
                    self.lambdalim = [6513, 6613]
                    self->redraw
                end
                'q': begin
                    case 1 of
                        self.keystate eq 1 or self.keystate eq 2: begin
                            self.keystate = 0
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Pixel mask modification cancelled.'
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'e': begin
                    case self.keystate of
                        0: begin
                            self.lambda1 = coords[0]
                            self.keystate = 1
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value="Press 'e' again to exclude wavelength region from pixel mask."
                        end
                        1: begin
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Updating database ...'
                            lambda1 = self.lambda1 < coords[0]
                            lambda2 = self.lambda1 > coords[0]
                            self.keystate = 0
                            scienceall = *self.science                            
                            w = where(scienceall[self.i].lambda/(1d + scienceall[self.i].z) gt lambda1 and scienceall[self.i].lambda/(1d + scienceall[self.i].z) lt lambda2, c)
                            if c eq 0 then begin
                                widget_control, widget_info(self.base, find_by_uname='status'), set_value='This wavelength region is invalid.'
                            endif else begin                            
                                scienceall[self.i].fitmask[w] = 0
                                ptr_free, self.science
                                self.science = ptr_new(scienceall)
                                self->redraw
                            endelse
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'i': begin
                    case self.keystate of
                        0: begin
                            self.lambda1 = coords[0]
                            self.keystate = 2
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value="Press 'i' again to include wavelength region in pixel mask."
                        end
                        2: begin
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Updating database ...'
                            lambda1 = self.lambda1 < coords[0]
                            lambda2 = self.lambda1 > coords[0]
                            self.keystate = 0
                            scienceall = *self.science                            
                            w = where(scienceall[self.i].lambda/(1d + scienceall[self.i].z) gt lambda1 and scienceall[self.i].lambda/(1d + scienceall[self.i].z) lt lambda2, c)
                            if c eq 0 then begin
                                widget_control, widget_info(self.base, find_by_uname='status'), set_value='This wavelength region is invalid.'
                            endif else begin                            
                                scienceall[self.i].fitmask[w] = 1
                                ptr_free, self.science
                                self.science = ptr_new(scienceall)
                                self->redraw
                            endelse
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase

                end
                'o': self->output_objname
                else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
            endcase
        end
    endcase
end


; ============== DESTROY ===============
pro sps_fit_cleanup, ev    
    widget_control, ev, get_uvalue=obj
    obj_destroy, obj
end


pro sps_fit::cleanup
    ptr_free, self.science
end


pro sps_fit::reprepare_all
    update_phot = 1

    curi = self.i
    for i=0,self.nspec-1 do begin
        self.i = i
        if ~update_phot and (*self.science)[self.i].good eq 0 then continue
        self->reprepare, /nostatusbar
    endfor
    scienceall = *self.science
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    self->writescience
    self->specres_mask, self.directory
    self.i = curi
    science = (*self.science)[self.i]
    self->statusbox, science=science
    self->redraw
end


pro sps_fit::default_cont
    scienceall = *self.science
    science = scienceall[self.i]
    science.contmask = 0
    self->continuum, science
    scienceall[self.i] = science
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    self->redraw
end


pro sps_fit::default_mask
    scienceall = *self.science
    science = scienceall[self.i]
    self->mask, science, nmc=0
    scienceall[self.i] = science
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    self->redraw
end

pro sps_fit::default_maskall
  curi = self.i
  for i=0,self.nspec-1 do begin
     self.i=i
     scienceall = *self.science
     science = scienceall[self.i]
     self->mask, science, nmc=0
     scienceall[self.i] = science
     ptr_free, self.science
     self.science = ptr_new(scienceall)
  endfor
  self.i = curi
  self->redraw
end

pro sps_fit::default_goodspec
    scienceall = *self.science
    curi = self.i
    for i=0,self.nspec-1 do begin
       self.i=i
       science = scienceall[self.i]
       sn = science.sn
       z  = science.z
       if sn gt 6. and z gt 0.05 then science.good = 1 else science.good = 0
       scienceall[self.i] = science
    endfor
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    self->writescience
    self.i = curi
    science = scienceall[self.i]
    self->statusbox, science=science
    self->redraw

end



pro sps_fit::reprepare, nostatusbar=nostatusbar
    update_phot = 1
    update_else = 1

    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Repreparing ...'
    scienceall = *self.science
    science = scienceall[self.i]
    contmask = science.contmask

    n = n_elements(science.lambda)
    wcont = where(contmask[3:n-4] eq 1)+3
    wcont = wcont[where(finite(science.telldiv[wcont]) and finite(science.continuum[wcont]) and science.continuum[wcont] ne 0)]
    dev = abs((science.telldiv[wcont] - science.continuum[wcont]) / science.continuum[wcont])
    avgdev = mean(dev)
    w = where(dev lt 3.0*avgdev, c)
    if c gt 0 then science.sn = 1.0/mean(dev[w])

    if update_phot eq 1 then begin
        common mask_in, mask_in
        case 1 of
            strmid(mask_in, 0, 4) eq '0024' or strmid(mask_in, 0, 3) eq 'rse': begin
                photfile = '/scr2/nichal/keck/deimos/Cl0024MS0451/Cl0024master.v7.fits.gz'
                vdispfile = '/scr2/nichal/keck/deimos/Cl0024MS0451/Cl0024_sigmas.txt'
                readcol, vdispfile, vd_objname, vd_vdisp, format='A,D', comment='#', /silent
            end
            strmid(mask_in, 0, 4) eq '0451': begin
                photfile = '/scr2/nichal/keck/deimos/Cl0024MS0451/MS0451master.v14.fits.gz'
                vdispfile ='/scr2/nichal/keck/deimos/Cl0024MS0451/MS0451_sigmas.txt'
                readcol, vdispfile, vd_objname, vd_vdisp, vd_vdisperr, format='A,D,D', comment='#', /silent
            end
            strmid(mask_in, 0, 6) eq 'Cl1604': begin
                photfile = getenv('UCI')+'Cl1604/catalogs/Cl1604.fits.gz'
            end
            else: message, mask_in+' is not associated with a cluster that I know.'
        endcase
        phot = mrdfits(photfile, 1, /silent)
        phot = phot[where(phot.dec gt -90 and phot.dec lt 90)]

        if strlowcase(strmid(mask_in, 0, 6)) eq 'cl1604' then begin
            match, strlowcase(strtrim(phot.mask, 2))+'_'+strtrim(phot.slit, 2), strlowcase(strtrim(science.mask, 2))+'_'+strtrim(science.slit, 2), w1, w2
        endif else begin
            spherematch, phot.ra, phot.dec, science.ra, science.dec, 1.0/3600., w1, w2
        endelse
        if w1[0] eq -1 then begin
        endif else begin
            science.z = phot[w1].z
            science.zquality = phot[w1].zquality
            if strlowcase(strmid(mask_in, 0, 6)) eq 'cl1604' then begin
                science.zspec = 0.0
                science.ra = phot[w1].ra
                science.dec = phot[w1].dec
                science.r = phot[w1].rmag
                science.i = phot[w1].imag
                science.f606w = phot[w1].f606wmag
                science.f814w = phot[w1].f814wmag
            endif else begin
                science.zspec = science.z
                science.zsource = phot[w1].zsource
                science.b = phot[w1].b_auto
                science.v = phot[w1].v_auto
                science.r = phot[w1].r_auto
                science.i = phot[w1].i_auto
                science.j = phot[w1].j_auto
                science.k = phot[w1].k_auto
                science.f814w = phot[w1].f814w_auto
                science.berr = phot[w1].b_auto_err
                science.verr = phot[w1].v_auto_err
                science.rerr = phot[w1].r_auto_err
                science.ierr = phot[w1].i_auto_err
                science.jerr = phot[w1].j_auto_err
                science.kerr = phot[w1].k_auto_err
                science.f814werr = phot[w1].f814w_auto_err
            endelse

            science.vdisp_smm = -999d
            science.vdisperr_smm = -999d
            if (size(vd_objname))[0] gt 0 then begin
                match, strtrim(vd_objname, 2), strtrim(science.objname, 2), w1, w2, count=cmatch
                if cmatch gt 0 then begin
                    science[w2].vdisp_smm = vd_vdisp[w1]
                    if strmid(mask_in, 0, 4) eq '0451' then begin
                        science[w2].vdisperr_smm = vd_vdisperr[w1]
                    endif else begin
                        science[w2].vdisperr_smm = -999d
                    endelse
                endif
            endif
        endelse
        case 1 of
            science.b gt 10.0 and science.b lt 50.0 and science.v gt 10.0 and science.v lt 50.0: science.phot_color = 'BV'
            science.b gt 10.0 and science.b lt 50.0 and science.r gt 10.0 and science.r lt 50.0: science.phot_color = 'BR'
            science.v gt 10.0 and science.v lt 50.0 and science.i gt 10.0 and science.i lt 50.0: science.phot_color = 'VI'
            science.v gt 10.0 and science.v lt 50.0 and science.k gt 10.0 and science.k lt 50.0: science.phot_color = 'VK'
            science.v gt 10.0 and science.v lt 50.0 and science.j gt 10.0 and science.j lt 50.0: science.phot_color = 'VJ'
            science.j gt 10.0 and science.j lt 50.0 and science.k gt 10.0 and science.k lt 50.0: science.phot_color = 'JK'
            science.f606w gt 10.0 and science.f606w lt 50.0 and science.f814w gt 10.0 and science.f814w lt 50.0: science.phot_color = 'F606WF814W'
            else: science.phot_color = 'BV'
        endcase
    endif
    
    if update_else eq 1 then begin
        ;self->skytweak, science
        ;science.skylinemask = -1
        ;science.contmask = 0
        self->specres, science, /goodoverride
        self->telluric, science
        self->continuum, science
        self->mask, science, /nomask
        science.spscont = 1.0
        self->indices, science, /noredraw
        ;self->fit, science, /noredraw, nostatusbar=nostatusbar
    endif

    self->statusbox, science=science
    scienceall[self.i] = science
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    self->redraw
end


; ============== SPECRES ==============
pro fitskyspec_gauss, x, a, f, pder
    u = (x - a[0])/a[1]
    f = a[2]*exp(-0.5*u^2)
    pder = [[(a[0]-x)/a[1]*f], [(x-a[0])^2./a[1]^3.*f], [f/a[2]]]
    return
end


pro sps_fit::skytweak, science
    blue = mrdfits(strtrim(science.spec1dfile, 2), 1, headblue, /silent)
    red = mrdfits(strtrim(science.spec1dfile, 2), 2, headred, /silent)
    lambdablue = blue.lambda
    lambdared = red.lambda
    path = file_dirname(science.spec1dfile)
    lambdabluenew = applytweaks(lambdablue, headblue, path)
    lambdarednew = applytweaks(lambdared, headred, path)
    lambdanew = [lambdabluenew, lambdarednew]
    science.lambda = lambdanew
end


function sps_fit::fitskylines, science
    lambda = science.lambda
    skyspec = science.skyspec
    ;skyspec = 1./science.ivar
    w = where(~finite(skyspec), c)
    if c gt 0 then skyspec[w] = 0
    skyspec = skyspec/max(skyspec)
    medskyspec = median(skyspec)
    ;plot, lambda, skyspec
    ;oplot, minmax(lambda), 2.0*replicate(medskyspec, 2), linestyle=1
    ;stop

    deriv1skyspec = deriv(lambda, skyspec)
    deriv2skyspec = deriv(lambda, deriv1skyspec)
    ;plot, lambda, deriv1skyspec
    ;oplot, minmax(lambda), replicate(0.2, 2), linestyle=1
    ;stop
    ;plot, lambda, deriv2skyspec, yrange=[-0.1, 0.1]
    ;oplot, minmax(lambda), replicate(-0.01, 2), linestyle=1
    ;stop
    thresh = 1.5
    nlines = 1000
    while nlines gt 200 do begin
        w = where(abs(deriv1skyspec) lt 0.2 and deriv2skyspec lt -0.01 and skyspec gt thresh*medskyspec)
        w = [w, n_elements(lambda)+1]
        wstep = round(-1*ts_diff(w, 1))
        linestart = where(wstep gt 1, nlines)
        if nlines lt 5 then begin
            message, 'Fewer than 5 sky lines in this spectrum.', /info
            return, [[-1], [-1], [-1]]
        endif
        linepix = round(-1*ts_diff(linestart, 1))
        nlines -= 1
        thresh *= 2
    endwhile
    linepix = linepix[0:nlines-1]
    wlocmax = lindgen(nlines)
    sigma = dblarr(nlines)
    sigmaerr = dblarr(nlines)
    linelambda = dblarr(nlines)
    nskyspec = n_elements(skyspec)
    wgoodline = bytarr(nlines)+1
    for i=0,nlines-1 do begin        
        if linepix[i] eq 1 then wlocmax[i] = w[linestart[i]] else begin
            junk = min(abs(deriv1skyspec[w[linestart[i]]:w[linestart[i]]+linepix[i]-1]), wmin)
            wlocmax[i] = w[linestart[i]]+wmin
        endelse
        if wlocmax[i]-10 lt 0 or wlocmax[i]+10 ge nskyspec then begin
            wgoodline[i] = 0
            continue
        endif
        skyspecmax = max(skyspec[wlocmax[i]-10:wlocmax[i]+10], wlocmaxtemp)
        wlocmax[i] += wlocmaxtemp-10
        wfit = lindgen(20)+wlocmax[i]-10
        wfit = wfit[where(wfit ge 0 and wfit lt n_elements(lambda), nfit)]
        lambdafit = lambda[wfit]
        skyspecfit = skyspec[wfit]

        ;lambdarange = max(lambda[wfit])-min(lambda[wfit])
        ;lambdastep = 0.03
        ;lambdafit = dindgen(round(lambdarange/lambdastep))*lambdastep+min(lambda[wfit])
        ;skyspecfit = interpol(skyspec[wfit], lambda[wfit], lambdafit, /quadratic)

        skyspecmax = max(skyspecfit, wmax)

        ;a = [lambdafit[wmax], 1.3/2.35, skyspecmax]
        ;yfit = curvefit(lambdafit, skyspecfit, lindgen(n_elements(lambdafit))+1.0, a, aerr, function_name='fitskyspec_gauss', chisq=chisq, status=status, /double, tol=1d-12)
        ;a = [skyspecmax, lambdafit[wmax], 1.3/2.35]

        guess = [skyspecmax, lambdafit[wmax], 0, medskyspec]
        yfit = gaussfit(lambdafit, skyspecfit, a, estimates=guess, sigma=aerr, nterms=4, chisq=chisq)
        sigma[i] = abs(a[2])
        sigmaerr[i] = aerr[2]
        linelambda[i] = a[1]

        ;print, sigma[i], sigmaerr[i]/sigma[i], chisq
        ;plot, lambdafit, skyspecfit
        ;oplot, lambdafit, yfit, color=fsc_color('red')
        ;wait, 0.5

        if chisq gt 1d-3 or sigmaerr[i] gt 0.8 then wgoodline[i] = 0
    endfor
    wgood = where(sigma gt 0 and wgoodline eq 1, c)
    if c gt 0 then begin
        linelambda = linelambda[wgood]
        sigma = sigma[wgood]
        sigmaerr = sigmaerr[wgood]
        ;ploterr, linelambda, sigma, sigmaerr, psym=1
        
        ;plot, lambda, skyspec
        ;oplot, lambda[wlocmax], skyspec[wlocmax], psym=1, color=fsc_color('red')

        return, [[linelambda], [sigma], [sigmaerr]]
    endif else return, [[-1], [-1], [-1]]
end


pro sps_fit::specres, science, qf=qf, goodoverride=goodoverride
    lambda = science.lambda
    wsky = where(science.skylinemask ne -1, cw)
    if cw eq 0 then begin
        fit = self->fitskylines(science)
        n = (size(fit))[1]
        if n gt 200 then message, 'Too many sky lines!'
        if n gt 0 then science.skyfit[0:n-1,*] = fit
        if n le 3 then begin
            science.dlam = replicate(1.37/2.35, n_elements(lambda))
            science.skylinemask = lonarr(n_elements(science.skylinemask))-1
            science.goodsky = 0
            message, 'Unstable sky line fit.  Using FWHM = 1.37 A', /info
            return
        endif

        w = where(2.35*fit[*,1] gt 0.8 and 2.35*fit[*,1] lt 7.0, cwprev)
        if cwprev lt 3 then begin
            science.dlam = replicate(1.37/2.35, n_elements(lambda))
            science.skylinemask = lonarr(n_elements(science.skylinemask))-1
            science.goodsky = 0
            message, 'Unusuable arc lines.  Using FWHM = 1.37 A', /info
            return
        endif

        ;quadratic fit
        qf = poly_fit(fit[w,0]/1000.0 - 7.8, 2.35*fit[w,1], 2, measure_errors=2.35*fit[w,2], chisq=chisq, /double, yfit=yfit)

        for j=0,4 do begin
            wnow = where(abs(2.35*fit[w,1] - yfit) lt 2.*2.35*fit[w,2], cw)
            if cw eq cwprev then break
            cwprev = cw
            if cw lt 3 then begin
                science.goodsky = 0
                message, 'The spectral resolution fit is very poor.', /info
                break
            endif
            w = w[wnow]
            qf = poly_fit(fit[w,0]/1000.0 - 7.8, 2.35*fit[w,1], 2, measure_errors=2.35*fit[w,2], chisq=chisq, /double, yfit=yfit)
        endfor
        n = (size(fit))[1]
        science.skylinemask = 0
        science.skylinemask[w] = 1
        if n lt 200 then science.skylinemask[n:n_elements(science.skylinemask)-1] = -1
    endif else begin
        wsky = where(science.skylinemask eq 1, csky)
        if csky lt 3 then begin
            science.dlam = replicate(1.37/2.35, n_elements(lambda))
            science.goodsky = 0
            message, 'Too few arc lines for new fit.', /info
            return
        endif
        fit = science.skyfit[wsky,*]
        qf = poly_fit(fit[wsky,0]/1000.0 - 7.8, 2.35*fit[wsky,1], 2, measure_errors=2.35*fit[wsky,2], chisq=chisq, /double, yfit=yfit)        
    endelse

    l = lambda / 1000. - 7.8
    dlam = poly(l, qf)
    dlam /= 2.35
    science.dlam = dlam
    if ~keyword_set(goodoverride) then science.goodsky = 1
end


pro sps_fit::add_skyline, lambda
    science = (*self.science)[self.i]
    wsky = where(science.skylinemask ne -1)
    fit = science.skyfit[wsky,*]
    if fit[0,0] eq -1 then begin
        widget_control, widget_info(self.base, find_by_uname='status'), set_value='There are no available sky lines.'
        return
    endif
    w = where(fit[*,0] ge self.lambdalim[0] * (1d + science.zspec) and fit[*,0] le self.lambdalim[1] * (1d + science.zspec), c)
    if c eq 0 then begin
        widget_control, widget_info(self.base, find_by_uname='status'), set_value='There are no available sky lines in this range of wavelengths.'
        return
    endif
    junk = min(abs(fit[w,0] - lambda), wmin)
    i = w[wmin]
    skylinemask = where(science.skylinemask[wsky] eq 1)
    if contains(skylinemask, i) then begin
        widget_control, widget_info(self.base, find_by_uname='status'), set_value='This sky line at '+string(round(fit[i,0]), format='(I4)')+' A is already included in the fit.'
        return
    endif else begin
        science.skylinemask[i] = 1
        self->specres, science
        scienceall = *self.science
        ptr_free, self.science
        scienceall[self.i] = science
        self.science = ptr_new(scienceall)
        self->redraw
    endelse
end


pro sps_fit::delete_skyline, lambda
    science = (*self.science)[self.i]
    wsky = where(science.skylinemask ne -1)
    fit = science.skyfit
    if fit[0,0] eq -1 then begin
        widget_control, widget_info(self.base, find_by_uname='status'), set_value='There are no available sky lines.'
        return
    endif
    w = where(fit[*,0] ge self.lambdalim[0] * (1d + science.zspec) and fit[*,0] le self.lambdalim[1] * (1d + science.zspec), c)
    if c eq 0 then begin
        widget_control, widget_info(self.base, find_by_uname='status'), set_value='There are no available sky lines in this range of wavelengths.'
        return
    endif
    junk = min(abs(fit[w,0] - lambda), wmin)
    i = w[wmin]
    skylinemask = where(science.skylinemask[wsky] eq 1)
    if contains(skylinemask, i) then begin
        w = where(skylinemask ne i, c)
        if c gt 0 then begin
            science.skylinemask[i] = 0
            self->specres, science
            scienceall = *self.science
            ptr_free, self.science
            scienceall[self.i] = science
            self.science = ptr_new(scienceall)
            self->redraw
        endif else begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='This sky line at '+string(round(fit[i,0]), format='(I4)')+' A is the only line in the fit.'
            return
        endelse
    endif else begin
        widget_control, widget_info(self.base, find_by_uname='status'), set_value='This sky line at '+string(round(fit[i,0]), format='(I4)')+' A is already excluded from the fit.'
        return
    endelse
end


pro sps_fit::specres_mask, directory
    sps_fit = mrdfits(directory+'/sps_fit.fits.gz', 1, /silent)
    nspec = n_elements(sps_fit)
    wuse = lonarr(nspec)
    d = 0
    for i=0,nspec-1 do begin
        w = where(sps_fit[i].skylinemask eq 1, nsky)
        if nsky gt 20 then begin
            lambda = d eq 0 ? reform(sps_fit[i].skyfit[w,0]) : [lambda, reform(sps_fit[i].skyfit[w,0])]
            res = d eq 0 ? reform(sps_fit[i].skyfit[w,1]) : [res, reform(sps_fit[i].skyfit[w,1])]
            err = d eq 0 ? reform(sps_fit[i].skyfit[w,1]) : [err, reform(sps_fit[i].skyfit[w,2])]
            d = 1
        endif
    endfor
    if (size(res))[0] gt 0 then begin
        w = where(2.35*res gt 0.8 and 2.35*res lt 7.0 and 2.35*err lt 1.0)
        specres_poly = poly_fit(lambda[w]/1000 - 7.8, 2.35*res[w], 2, measure_errors=2.35*err[w], /double)
        specres_lambda = dindgen(8192)*(9100-6300)/8191 + 6300
        ;ploterror, lambda[w], 2.35*res[w], 2.35*err[w], /nohat, psym=1
        ;oplot, specres_lambda, poly(specres_lambda/1000 - 7.8, specres_poly), color=fsc_color('red')
    endif else specres_poly = [2.35, 0.0]
    save, specres_poly, filename=directory+'/specres_poly.sav'
end


; ============= CONTINUUM =============
pro sps_fit::continuum, science
    fft = 0
    spline = 1
    poly = 0
    usesmooth = 0

    contmask = science.contmask
    w = where(contmask ne 0, c)
    lambda = science.lambda / (1d + science.zspec)
    n = n_elements(lambda)
    nhalf = round(double(n) / 2.)
    if c eq 0 then begin
        linestart = *self.linestart
        lineend = *self.lineend
        contmask = bytarr(n_elements(lambda))+1
        for i=0,n_elements(linestart)-1 do begin
            w = where(lambda ge linestart[i] and lambda le lineend[i], c)
            if c gt 0 then contmask[w] = 0
        endfor
        tellmask = bytarr(n_elements(lambda))
        tellstart = [6864., 7591., 8938.]
        tellend = [6935., 7694., 100000.]
        for i=0,n_elements(tellstart)-1 do begin
            w = where(science.lambda ge tellstart[i] and science.lambda le tellend[i], c)
            if c gt 0 then begin
                contmask[w] = 0
                tellmask[w] = 1
            endif
        endfor
        contmask[0:2] = 0
        contmask[nhalf-3:nhalf+2] = 0
        contmask[n-3:n-1] = 0
     endif

    satbands = [6870, 7650]
    niter = 5

    if n lt 8193 then begin
        wwhole = lindgen(nhalf)
        ccd2 = 2
    endif else begin
        wwhole = lindgen(n)
        ccd2 = 1
    endelse
    spec = science.telldiv
    for ccd=1,ccd2 do begin
        contmask[wwhole[0:3]] = 0
        contmask[wwhole[nhalf-4:nhalf-1]] = 0
        won = where(contmask[wwhole] eq 1, complement=woff, con) + wwhole[0]
        woff += wwhole[0]
        
        case 1 of
            spline: begin
                bkpt = slatec_splinefit(lambda[won], spec[won], coeff, invvar=science.telldivivar[won], bkspace=150, upper=3, lower=3, /silent)
                if bkpt[0] eq -1 then return
                cont = slatec_bvalu(lambda[wwhole], bkpt, coeff)
            end
            poly: begin
                degree = 12
                norm = median(spec[won])
                a = [norm, replicate(0.0, degree-1)]
                p = lmfit(lambda[won], spec[won], a, measure_errors=(science.telldivivar[won])^(-0.5), /double, function_name='legendre_poly')
                cont = legendre_poly(lambda[wwhole], a, /noderiv)
            end
            fft: begin
                wfft = wwhole[10:nhalf-11]
                nfft = n_elements(wfft)
                nkeep = 10
                ft = fft(spec[wfft], -1, /double)

                ft[nkeep:nfft-1] = 0
                cont = real_part(fft(ft, 1, /double))
                plot, cont
                stop
            end
            usesmooth: begin
                ww = wwhole[3:nhalf-4]
                wcont = where(contmask[ww] eq 1, ccont) + ww[0]
                for i=1,5 do begin
                    wcontold = wcont
                    if ccont lt 50 then begin
                        science.continuum = replicate(-999, n_elements(science.continuum))
                        return
                    endif
                    cont = smooth_gauss_wrapper(lambda[wcont], spec[wcont], lambda[wwhole], 30.0, ivar1=science.telldivivar[wcont])
                    wcont = where(abs(spec[ww]-cont[ww-wwhole[0]]) lt (science.telldivivar[ww])^(-0.5) and contmask[ww] eq 1, ccont) + ww[0]

                    if array_equal(wcont, wcontold) then break
                endfor
            end
        endcase

        if ccd eq 1 then contb = cont
        if ccd eq 2 then contr = cont
        wwhole += nhalf
    endfor
    if n lt 8193 then cont = [contb, contr] else cont = contb

    science.contmask = contmask
    science.continuum = cont
    science.contdiv = science.telldiv / cont
    science.contdivivar = science.telldivivar * cont^2.

    wcont = where(contmask[3:n-4] eq 1)+3
    wcont = wcont[where(finite(science.telldiv[wcont]) and finite(science.continuum[wcont]) and science.continuum[wcont] ne 0)]
    dev = abs((science.telldiv[wcont] - science.continuum[wcont]) / science.continuum[wcont])
    avgdev = mean(dev)
    w = where(dev lt 3.0*avgdev, c)
    if c gt 0 then science.sn = 1.0/mean(dev[w])
end


; ============= TELLURIC ==============
pro sps_fit::telluric, science
    specresfile = self.directory+'specres_poly.sav'
    globalres = 0
    if file_test(specresfile) then begin
        restore, specresfile
        dlam = poly(science.lambda/1000 - 7.8, specres_poly) / 2.35
        globalres = 1
    endif

    aratio = science.airmass/((*self.tell).airmass)
    telllambda = (*self.tell)[0].lambda
    tellspec = ((*self.tell)[0].spec)^aratio
    tellivar = ((*self.tell)[0].ivar)*(((*self.tell)[0].spec)/(aratio*tellspec))^2.
    ivarmissing = 10d10
    w = where(tellivar ge 1d8, c)
    if c gt 0 then tellivar[w] = ivarmissing
    f = where(finite(tellspec) and tellspec gt 0 and finite(tellivar) and tellivar gt 0 and tellivar lt 1d8)
    telllambda = telllambda[f]
    tellspec = tellspec[f]
    tellivarfull = tellivar
    tellivar = tellivar[f]

    tellspecnew = interpolate(tellspec, findex(telllambda, science.lambda), missing=1.)
    tellivarnew = interpolate(tellivarfull, findex(telllambda, science.lambda), missing=ivarmissing)
    science.telldiv = science.spec / tellspecnew
    science.telldivivar = (science.ivar^(-1.) + (science.telldiv)^2.*tellivarnew^(-1.))^(-1.) * tellspecnew^2.
end


; =============== MASK ================
pro sps_fit::mask, science, nomask=nomask, zfind=zfind, nozfind=nozfind, nmc=nmc
    specresfile = self.directory+'specres_poly.sav'
    globalres = 0
    if file_test(specresfile) then begin
        restore, specresfile
        dlam = poly(science.lambda/1000 - 7.8, specres_poly) / 2.35
        globalres = 1
    endif

    if ~keyword_set(nomask) then begin
        linestart = *self.linestart
        lineend = *self.lineend
        linetype = *self.linetype
        wem = where(linetype eq 'e', cemlines)
        mask = bytarr(n_elements(science.lambda))+1
        w = where(science.ivar le 0 or ~finite(science.ivar) or science.lambda/(1d + science.zspec) lt 3650, cw)
        if cw gt 0 then mask[w] = 0
        ;for i=0,cemlines-1 do begin
        ;    w = where(science.lambda/(1d + science.zspec) ge linestart[wem[i]] and science.lambda/(1d + science.zspec) le lineend[wem[i]], c)
        ;    if c gt 0 then mask[w] = 0
        ;endfor
        tellstart = [6864., 7591., 8938.]
        tellend = [6935., 7694., 100000.]
        for i=0,n_elements(tellstart)-1 do begin
            w = where(science.lambda ge tellstart[i] and science.lambda le tellend[i], c)
            if c gt 0 then mask[w] = 0
        endfor
        w = where((science.contdiv-science.continuum) gt 3.*science.contdivivar^(-0.5), cw)
        if cw gt 0 then mask[w] = 0
        n = n_elements(science.lambda)
        nhalf = round(double(n) / 2.)
        mask[0:4] = 0
        mask[nhalf-4:nhalf+4] = 0
        mask[n-5:n-1] = 0
    endif else begin
        mask = science.fitmask
    endelse
    science.fitmask = mask
end


; =============== MASK Fe bands================
pro sps_fit::maskfe, science, nomask=nomask, zfind=zfind, nozfind=nozfind, nmc=nmc
    specresfile = self.directory+'specres_poly.sav'
    globalres = 0
    if file_test(specresfile) then begin
        restore, specresfile
        dlam = poly(science.lambda/1000 - 7.8, specres_poly) / 2.35
        globalres = 1
    endif

    if ~keyword_set(nomask) then begin
        linestart = *self.linestart
        lineend  = *self.lineend
        linetype = *self.linetype
        indstart = *self.indstart 
        indend   = *self.indend   
        indname  = *self.indname 
        mask = bytarr(n_elements(science.lambda))

        ;mask everywhere else = 0 except where the indices bands are
        use_indices = ['Fe4383','Fe4531','Fe4668','Fe5015','Fe5270','Fe5335','Fe5406','Fe5709','Fe5782']
        for i=0,n_elements(use_indices)-1 do begin
           indnow = where(indname eq use_indices(i),cindnow)
           if cindnow eq 0 then stop
           w = where(science.lambda/(1.+science.zspec) gt indstart(indnow[0]) and science.lambda/(1.+science.zspec) lt indend(indnow[0]),cw)
           if cw gt 0 then mask[w]=1
        endfor

        w = where(science.ivar le 0 or ~finite(science.ivar) or science.lambda/(1d + science.zspec) lt 3650 or science.lambda/(1.+science.zspec) gt 7400, cw)
        if cw gt 0 then mask[w] = 0 ;bad points
        w = where(science.fitmask eq 0, cw)
        if cw gt 0 then mask[w] = 0 ;bad points
        tellstart = [6864., 7591., 8938.]
        tellend = [6935., 7694., 100000.]
        for i=0,n_elements(tellstart)-1 do begin
            w = where(science.lambda ge tellstart[i] and science.lambda le tellend[i], c)
            if c gt 0 then mask[w] = 0
        endfor
        w = where((science.contdiv-science.continuum) gt 3.*science.contdivivar^(-0.5), cw)
        if cw gt 0 then mask[w] = 0
        n = n_elements(science.lambda)
        mask[0:4] = 0
        mask[n-5:n-1] = 0
    endif else begin
        mask = science.fitfemask
    endelse
    science.fitfemask = mask
 end

; =============== MASK Mg bands================
pro sps_fit::maskmg, science, nomask=nomask, zfind=zfind, nozfind=nozfind, nmc=nmc
    specresfile = self.directory+'specres_poly.sav'
    globalres = 0
    if file_test(specresfile) then begin
        restore, specresfile
        dlam = poly(science.lambda/1000 - 7.8, specres_poly) / 2.35
        globalres = 1
    endif

    if ~keyword_set(nomask) then begin
        linestart = *self.linestart
        lineend  = *self.lineend
        linetype = *self.linetype
        indstart = *self.indstart 
        indend   = *self.indend   
        indname  = *self.indname 
        mask = bytarr(n_elements(science.lambda))

        ;mask everywhere else = 0 except where the indices bands are
        use_indices = ['Mg_b','Mg_2','Mg_1']
        for i=0,n_elements(use_indices)-1 do begin
           indnow = where(indname eq use_indices(i),cindnow)
           if cindnow eq 0 then stop
           w = where(science.lambda/(1.+science.zspec) gt indstart(indnow[0]) and science.lambda/(1.+science.zspec) lt indend(indnow[0]),cw)
           if cw gt 0 then mask[w]=1
        endfor

        w = where(science.ivar le 0 or ~finite(science.ivar) or science.lambda/(1d + science.zspec) lt 3650 or science.lambda/(1.+science.zspec) gt 7400, cw)
        if cw gt 0 then mask[w] = 0 ;bad points
        w = where(science.fitmask eq 0, cw)
        if cw gt 0 then mask[w] = 0 ;bad points

        tellstart = [6864., 7591., 8938.]
        tellend = [6935., 7694., 100000.]
        for i=0,n_elements(tellstart)-1 do begin
            w = where(science.lambda ge tellstart[i] and science.lambda le tellend[i], c)
            if c gt 0 then mask[w] = 0
        endfor
        w = where((science.contdiv-science.continuum) gt 3.*science.contdivivar^(-0.5), cw)
        if cw gt 0 then mask[w] = 0
        n = n_elements(science.lambda)
        mask[0:4] = 0
        mask[n-5:n-1] = 0
    endif else begin
        mask = science.fitmgmask
    endelse
    science.fitmgmask = mask
end

; ============== REDRAW ===============
pro sps_fit::redraw
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Redrawing ...'
    self->default_range, /update
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    widget_control, widget_info(self.base, find_by_uname='spec'), get_value=index
    wset, index
    science = (*self.science)[self.i]

    case mode of
        -1: begin                ;smoothed
            specsmooth = smooth(science.spec, 100)
            plot, science.lambda, specsmooth, xrange=self.lambdalim * (1d + science.zspec), yrange=self.ylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6100-pixel boxcar smoothed flux (e!E-!N/hr)!3'
            yoff1 = 0.01*(self.ylim[1] - self.ylim[0])
            yoff2 = 0.05*(self.ylim[1] - self.ylim[0])
            n = n_elements(*self.tellstart)
            for i=0,n-1 do begin
                oplot, [(*self.tellstart)[i], (*self.tellend)[i]], 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(*self.tellthick)[i]
            endfor
            n = n_elements(*self.linewaves)
            for i=0,n-1 do begin
                if (*self.linewaves)[i] * (1d + science.zspec) le !X.CRANGE[0] or (*self.linewaves)[i] * (1d + science.zspec) ge !X.CRANGE[1] then continue
                oplot, [(*self.linewaves)[i], (*self.linewaves)[i]] * (1d + science.zspec), [0.06*!Y.CRANGE[0]+0.94*!Y.CRANGE[1], 0.02*!Y.CRANGE[0]+0.98*!Y.CRANGE[1]], color=fsc_color((*self.linecolors)[i])
                xyouts, ((*self.linewaves)[i] * (1d + science.zspec))+0.002*(!X.CRANGE[1]-!X.CRANGE[0]), 0.07*!Y.CRANGE[0]+0.93*!Y.CRANGE[1], (*self.linenames)[i], orientation=90, alignment=1, color=fsc_color((*self.linecolors)[i])
            endfor
        end

        0: begin        ;telluric cross-correlation
            tell = (*self.tell)[0]
            w = where(tell.spec gt 0)
            cont = interpolate(science.continuum, findex(science.lambda/(1d + science.zspec), tell.lambda[w]))
            plot, science.lambda/(1d + science.zspec), science.spec, xrange=self.lambdalim * (1d + science.zspec), yrange=self.ylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6flux (e!E-!N/hr)!3'
            oplot, tell.lambda[w], tell.spec[w]*cont, color=fsc_color('green')
            n = n_elements(*self.tellstart)
            for i=0,n-1 do begin
                oplot, [(*self.tellstart)[i], (*self.tellend)[i]], 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(*self.tellthick)[i]
            endfor
        end

        1: begin                ;continuum fit
            t = round(-1*ts_diff(science.contmask, 1))
            wstart = where(t eq 1, cstart)+1
            wend = where(t eq -1, cend)
            if science.contmask[0] eq 1 then begin
                if cstart eq 0 then begin
                    wstart = 0
                endif else begin
                    wstart = [0, wstart]
                endelse
                cstart += 1
            endif
            if science.contmask[n_elements(t)-1] eq 1 then begin
                if cend eq 0 then begin
                    wend = n_elements(t)-1
                endif else begin
                    wend = [wend, n_elements(t)-1]
                endelse
                cend += 1
            endif
            if cstart ne cend then message, 'There are a different number of starting and ending continuum wavelengths.'
            if cstart eq 0 or cend eq 0 then message, 'There are no continuum regions.'

            plot, science.lambda, science.telldiv/median(science.telldiv), xrange=self.lambdalim * (1d + science.zspec), yrange=self.ylim, xstyle=5, ystyle=5, background=fsc_color('white'), color=fsc_color('black'), /nodata
            for i=0,cstart-1 do begin
                x = ([science.lambda[wstart[i]], science.lambda[wstart[i]], science.lambda[wend[i]], science.lambda[wend[i]]] > (self.lambdalim[0] * (1d + science.zspec))) < (self.lambdalim[1] * (1d + science.zspec))
                y = [self.ylim[0], self.ylim[1], self.ylim[1], self.ylim[0]]
                polyfill, x, y, color=fsc_color('light cyan')
            endfor
            oplot, science.lambda, science.telldiv/median(science.telldiv), color=fsc_color('black')
            oplot, science.lambda, science.continuum/median(science.telldiv), color=fsc_color('green')
            if science.feh lt 3. then begin
               normfactor =median(science.spsspecfull)
               oplot, science.lambda,science.spsspecfull/normfactor,color=fsc_color('red')            
               oplot, science.lambda,science.spscontfull/normfactor,color=fsc_color('orange')            
            endif
            if science.fe lt 3. then begin
               normfactor = median(science.spsspecfullfe)
               oplot, science.lambda,science.spsspecfullfe/normfactor,color=fsc_color('darkgreen')            
               oplot, science.lambda,science.spscontfullfe/normfactor,color=fsc_color('darkgreen')            
            endif
            if science.mg lt 3. then begin
               normfactor = median(science.spsspecfullmg)
               oplot, science.lambda,science.spsspecfullmg/normfactor,color=fsc_color('blue')            
               oplot, science.lambda,science.spscontfullmg/normfactor,color=fsc_color('blue')            
            endif

            plot, science.lambda, science.telldiv, xrange=self.lambdalim * (1d + science.zspec), yrange=self.ylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6flux (e!E-!N/hr)!3', /nodata, /noerase
            n = n_elements(*self.tellstart)
            for i=0,n-1 do begin
                oplot, [(*self.tellstart)[i], (*self.tellend)[i]], 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(*self.tellthick)[i]
            endfor
            n = n_elements(*self.linewaves)
            for i=0,n-1 do begin
                if (*self.linewaves)[i] * (1d + science.zspec) le !X.CRANGE[0] or (*self.linewaves)[i] * (1d + science.zspec) ge !X.CRANGE[1] then continue
                oplot, [(*self.linewaves)[i], (*self.linewaves)[i]] * (1d + science.zspec), [0.06*!Y.CRANGE[0]+0.94*!Y.CRANGE[1], 0.02*!Y.CRANGE[0]+0.98*!Y.CRANGE[1]], color=fsc_color((*self.linecolors)[i])
                xyouts, ((*self.linewaves)[i] * (1d + science.zspec))+0.002*(!X.CRANGE[1]-!X.CRANGE[0]), 0.07*!Y.CRANGE[0]+0.93*!Y.CRANGE[1], (*self.linenames)[i], orientation=90, alignment=1, color=fsc_color((*self.linecolors)[i])
            endfor
        end

        2: begin                ;continuum division
            plot, science.lambda, science.contdiv, xrange=self.lambdalim * (1d + science.zspec), yrange=self.divylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6flux (normalized)!3', /nodata
            oplot, self.lambdalim, [1.0, 1.0], color=fsc_color('orange')
            oplot, science.lambda, science.contdiv, color=fsc_color('black')
            n = n_elements(*self.tellstart)
            for i=0,n-1 do begin
                oplot, [(*self.tellstart)[i], (*self.tellend)[i]], 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(*self.tellthick)[i]
            endfor
        end

        3: begin        ;sky line fit
            wsky = where((*self.science)[self.i].skylinemask ne -1, csky)
            plot, science.lambda, 2.35*science.dlam, xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6sky line FWHM (!sA!r!u!9 %!6!n)!3', xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xrange=self.lambdalim * (1d + science.zspec), yrange=self.skyylim
            if csky gt 0 then begin
                fit = science.skyfit[wsky,*]
                if fit[0,0] ne -1 then begin
                    oploterror, fit[wsky,0], 2.35*fit[wsky,1], 2.35*fit[wsky,2], psym=1, color=fsc_color('black'), errcolor=fsc_color('black'), /nohat
                    wdel = where(science.skylinemask[wsky] eq 0, cdel)
                    if cdel gt 0 then plots, fit[wdel,0], 2.35*fit[wdel,1], color=fsc_color('red'), psym=7
                endif
            endif
            specresfile = self.directory+'specres_poly.sav'
            if file_test(specresfile) then begin
                restore, specresfile
                oplot, science.lambda, poly(science.lambda/1000 - 7.8, specres_poly), color=fsc_color('green')
            endif
        end

        4: begin        ;rest frame
            t = round(-1*ts_diff(science.fitmask, 1))
            wstart = where(t eq 1, cstart)+1
            wend = where(t eq -1, cend)
            if science.fitmask[0] eq 1 then begin
                if cstart eq 0 then begin
                    wstart = 0
                endif else begin
                    wstart = [0, wstart]
                endelse
                cstart += 1
            endif
            if science.fitmask[n_elements(t)-1] eq 1 then begin
                if cend eq 0 then begin
                    wend = n_elements(t)-1
                endif else begin
                    wend = [wend, n_elements(t)-1]
                endelse
                cend += 1
            endif
            if cstart ne cend then message, 'There are a different number of starting and ending fitmask wavelengths.'
            
            plot, science.lambda/(1d + science.zspec), science.contdiv, xrange=self.lambdalim, yrange=self.divylim, xstyle=5, ystyle=5, background=fsc_color('white'), color=fsc_color('black'), /nodata

            if cstart eq 0 or cend eq 0 then message, 'There are no fitmask regions.', /info else begin
                for i=0,cstart-1 do begin
                    x = ([science.lambda[wstart[i]], science.lambda[wstart[i]], science.lambda[wend[i]], science.lambda[wend[i]]]/(1d + science.zspec) > self.lambdalim[0]) < self.lambdalim[1]
                    y = [self.divylim[0], self.divylim[1], self.divylim[1], self.divylim[0]]
                    polyfill, x, y, color=fsc_color('light cyan')
                endfor
            endelse
            oplot, [-1d6, 1d6], [0.0, 0.0], color=fsc_color('pink')
            oplot, [-1d6, 1d6], [1.0, 1.0], color=fsc_color('pale green')
            if science.zfit ne 0 and finite(science.zfit) then oplot, science.lambda/(1d + science.zfit), science.contdiv, color=fsc_color('black') else oplot, science.lambda/(1d + science.zspec), science.contdiv, color=fsc_color('black')

            n = n_elements(*self.tellstart)
            for i=0,n-1 do begin
                oplot, [(*self.tellstart)[i], (*self.tellend)[i]] / (1d + science.zspec), 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(*self.tellthick)[i]
            endfor
            n = n_elements(*self.linewaves)
            for i=0,n-1 do begin
                if (*self.linewaves)[i] le !X.CRANGE[0] or (*self.linewaves)[i] ge !X.CRANGE[1] then continue
                oplot, [(*self.linewaves)[i], (*self.linewaves)[i]], [0.06*!Y.CRANGE[0]+0.94*!Y.CRANGE[1], 0.02*!Y.CRANGE[0]+0.98*!Y.CRANGE[1]], color=fsc_color((*self.linecolors)[i])
                xyouts, (*self.linewaves)[i]+0.002*(!X.CRANGE[1]-!X.CRANGE[0]), 0.07*!Y.CRANGE[0]+0.93*!Y.CRANGE[1], (*self.linenames)[i], orientation=90, alignment=1, color=fsc_color((*self.linecolors)[i])
            endfor
            plot, science.lambda/(1d + science.zspec), science.contdiv / science.spscont, xrange=self.lambdalim, yrange=self.divylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6rest wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6flux (normalized)!3', /nodata, /noerase
            if science.feh gt -10 and science.age gt 0.0 and total(science.spsspec gt 0.0) then begin
                oplot, science.lambda/(1d + science.zspec), science.spsspec, color=fsc_color('red')
             endif

            if science.fe gt -10 then begin
               oplot, science.lambda/(1d + science.zspec), science.spsspecfe, color=fsc_color('darkgreen')
               t = round(-1*ts_diff(science.fitfemask, 1))
               wstart = where(t eq 1, cstart)+1
               wend = where(t eq -1, cend)
               if science.fitfemask[0] eq 1 then begin
                  if cstart eq 0 then begin
                     wstart = 0
                  endif else begin
                     wstart = [0, wstart]
                  endelse
                  cstart += 1
               endif
               if science.fitfemask[n_elements(t)-1] eq 1 then begin
                  if cend eq 0 then begin
                     wend = n_elements(t)-1
                  endif else begin
                     wend = [wend, n_elements(t)-1]
                  endelse
                  cend += 1
               endif
               if cstart ne cend then message, 'There are a different number of starting and ending fitmask wavelengths.'
               for i=0,cstart-1 do begin
                  x = [science.lambda[wstart[i]],science.lambda[wend[i]]]/(1d + science.zspec)
                  y = [0.2, 0.2]
                  oplot, x, y, color=fsc_color('darkgreen'), thick=5
               endfor
            endif
            if science.mg gt -10 then begin
               oplot, science.lambda/(1d + science.zspec), science.spsspecmg, color=fsc_color('blue')
                  t = round(-1*ts_diff(science.fitmgmask, 1))
               wstart = where(t eq 1, cstart)+1
               wend = where(t eq -1, cend)
               if science.fitmgmask[0] eq 1 then begin
                  if cstart eq 0 then begin
                     wstart = 0
                  endif else begin
                     wstart = [0, wstart]
                  endelse
                  cstart += 1
               endif
               if science.fitmgmask[n_elements(t)-1] eq 1 then begin
                  if cend eq 0 then begin
                     wend = n_elements(t)-1
                  endif else begin
                     wend = [wend, n_elements(t)-1]
                  endelse
                  cend += 1
               endif
               if cstart ne cend then message, 'There are a different number of starting and ending fitmask wavelengths.'
               for i=0,cstart-1 do begin
                  x = [science.lambda[wstart[i]],science.lambda[wend[i]]]/(1d + science.zspec)
                  y = [0.15, 0.15]
                  oplot, x, y, color=fsc_color('blue'), thick=5
               endfor
            endif               
         end
    endcase        
    zl = mode lt 4 ? (*self.science)[self.i].z : 0.0
    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)

    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end


pro sps_fit::statusbox, science=science
    if ~keyword_set(science) then science = (*self.science)[self.i]
    
    phot_color = science.phot_color
    case phot_color of
        'VK': begin
            color = science.v-science.k
            mag = science.k
            colorerr = sqrt((science.verr)^2. + (science.kerr)^2.)
            magerr = science.kerr
            noerror = science.verr le 0 or science.verr ge 10 or science.kerr le 0 or science.kerr ge 10
            clrlbl = 'V-K = '
            maglbl = 'K = '
        end
        'VJ': begin
            color = science.v-science.j
            mag = science.j
            colorerr = sqrt((science.verr)^2. + (science.jerr)^2.)
            magerr = science.jerr
            noerror = science.verr le 0 or science.verr ge 10 or science.jerr le 0 or science.jerr ge 10
            clrlbl = 'V-J = '
            maglbl = 'J = '
        end
        'VI': begin
            color = science.v-science.i
            mag = science.i
            colorerr = sqrt((science.verr)^2. + (science.ierr)^2.)
            magerr = science.ierr
            noerror = science.verr le 0 or science.verr ge 10 or science.ierr le 0 or science.ierr ge 10
            clrlbl = 'V-I = '
            maglbl = 'I = '
        end
        'VR': begin
            color = science.v-science.r
            mag = science.r
            colorerr = sqrt((science.verr)^2. + (science.rerr)^2.)
            magerr = science.rerr
            noerror = science.verr le 0 or science.verr ge 10 or science.rerr le 0 or science.rerr ge 10
            clrlbl = 'V-R = '
            maglbl = 'R = '
        end
        'BV': begin
            color = science.b-science.v
            mag = science.v
            colorerr = sqrt((science.berr)^2. + (science.verr)^2.)
            magerr = science.verr
            noerror = science.berr le 0 or science.berr ge 10 or science.verr le 0 or science.verr ge 10
            clrlbl = 'B-V = '
            maglbl = 'V = '
        end
        'BR': begin
            color = science.b-science.r
            mag = science.r
            colorerr = sqrt((science.berr)^2. + (science.rerr)^2.)
            magerr = science.rerr
            noerror = science.berr le 0 or science.berr ge 10 or science.rerr le 0 or science.rerr ge 10
            clrlbl = 'B-R = '
            maglbl = 'R = '
        end
        'RI': begin
            color = science.r-science.i
            mag = science.i
            colorerr = sqrt((science.rerr)^2. + (science.ierr)^2.)
            magerr = science.ierr
            noerror = science.rerr le 0 or science.rerr ge 10 or science.ierr le 0 or science.ierr ge 10
            clrlbl = 'R-I = '
            maglbl = 'I = '
        end
        'JK': begin
            color = science.j-science.k
            mag = science.k
            colorerr = sqrt((science.jerr)^2. + (science.kerr)^2.)
            magerr = science.kerr
            noerror = science.jerr le 0 or science.jerr ge 10 or science.kerr le 0 or science.kerr ge 10
            clrlbl = 'J-K = '
            maglbl = 'K = '
        end
        'F606WF814W': begin
            color = science.f606w-science.f814w
            mag = science.f814w
            colorerr = sqrt((science.f606werr)^2. + (science.f814werr)^2.)
            magerr = science.f814werr
            noerror = science.f606werr le 0 or science.f606werr ge 10 or science.f814werr le 0 or science.f814werr ge 10
            clrlbl = 'F606W-F814W = '
            maglbl = 'F814W = '
        end
        else: begin
            color = -999.0
            mag = -999.0
            colorerr = 0.0
            magerr = 0.0
            noerror = 1
            clrlbl = 'V-I = '
            maglbl = 'I = '
        end
    endcase

    unknown = '???'
    widget_control, widget_info(self.base, find_by_uname='good'), set_value=[science.goodsky, science.good, science.goodfit]
    widget_control, widget_info(self.base, find_by_uname='curid'), set_value=strtrim(science.objname, 2)+' ('+strcompress(self.i+1, /rem)+' / '+strcompress(self.nspec, /rem)+')'
    widget_control, widget_info(self.base, find_by_uname='maglabel'), set_value=maglbl
    widget_control, widget_info(self.base, find_by_uname='collabel'), set_value=clrlbl
    widget_control, widget_info(self.base, find_by_uname='curmag'), set_value=mag gt 0 ? strcompress(string(mag, format='(D10.2)'), /rem)+(noerror ? '' : ' +/- '+strcompress(string(magerr, format='(D10.2)'), /rem)) : unknown
    widget_control, widget_info(self.base, find_by_uname='curcol'), set_value=color gt -10 ? strcompress(string(color, format='(D10.2)'), /rem)+(noerror ? '' : ' +/- '+strcompress(string(colorerr, format='(D10.2)'), /rem)) : unknown
    widget_control, widget_info(self.base, find_by_uname='curz'), set_value=strcompress(string(science.z, format='(D5.3)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='curzfit'), set_value=strcompress(string(science.zfit, format='(D5.3)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='curzquality'), set_value=strcompress(string(science.zquality, format='(D4.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='cursn'), set_value=science.sn gt 0 ? strcompress(string(science.sn, format='(D10.1)'), /rem) : unknown
    widget_control, widget_info(self.base, find_by_uname='curnloop'), set_value=science.nloop gt 0 ? strcompress(string(science.nloop, format='(D10.1)'), /rem) : unknown
    widget_control, widget_info(self.base, find_by_uname='curage'), set_value=science.age gt -100 ? strcompress(string(science.age, format='(D10.2)'), /rem)+(science.ageerr le 0 ? '' : ' +/- '+strcompress(string(science.ageerr, format='(D10.2)'), /rem))+' Gyr' : unknown
    widget_control, widget_info(self.base, find_by_uname='curmstar'), set_value=science.logmstar gt 0 ? strcompress(string(science.logmstar, format='(D10.2)'), /rem) : unknown
    widget_control, widget_info(self.base, find_by_uname='curfeh'), set_value=science.feh gt -100 ? strcompress(string(science.feh, format='(D10.2)'), /rem)+(science.feherr le 0 ? '' : ' +/- '+strcompress(string(science.feherr, format='(D10.2)'), /rem)) : unknown
    widget_control, widget_info(self.base, find_by_uname='curfe'), set_value=science.fe gt -100 ? strcompress(string(science.fe, format='(D10.2)'), /rem)+(science.feerr le 0 ? '' : ' +/- '+strcompress(string(science.feerr, format='(D10.2)'), /rem)) : unknown
    widget_control, widget_info(self.base, find_by_uname='curmg'), set_value=science.mg gt -100 ? strcompress(string(science.mg, format='(D10.2)'), /rem)+(science.mgerr le 0 ? '' : ' +/- '+strcompress(string(science.mgerr, format='(D10.2)'), /rem)) : unknown
    widget_control, widget_info(self.base, find_by_uname='curcah'), set_value=science.caherr gt 0 ? strcompress(string(science.cah, format='(D10.2)'), /rem)+(science.caherr le 0 ? '' : ' +/- '+strcompress(string(science.caherr, format='(D10.2)'), /rem)) : unknown
    widget_control, widget_info(self.base, find_by_uname='curgband'), set_value=science.gbanderr gt 0 ? strcompress(string(science.gband, format='(D10.2)'), /rem)+(science.gbanderr le 0 ? '' : ' +/- '+strcompress(string(science.gbanderr, format='(D10.2)'), /rem)) : unknown
    widget_control, widget_info(self.base, find_by_uname='curchisq'), set_value=science.chisq gt 0 ? strcompress(string(science.chisq, format='(D10.1)'), /rem) : unknown

    ;widget_control, widget_info(self.base, find_by_uname='curafe'), set_value=science.alphafe gt -100 ? strcompress(string(science.alphafe, format='(D10.2)'), /rem)+(science.alphafeerr le 0 ? '' : ' +/- '+strcompress(string(science.alphafeerr, format='(D10.2)'), /rem)) : unknown
    case 1 of
        science.vdisp gt 0: widget_control, widget_info(self.base, find_by_uname='curvdisp'), set_value=strcompress(string(science.vdisp, format='(D10.1)'), /rem)+(science.vdisperr le 0 ? '' : ' +/- '+strcompress(string(science.vdisperr, format='(D10.1)'), /rem))+' km/s'
        science.vdisp_smm gt 0: widget_control, widget_info(self.base, find_by_uname='curvdisp'), set_value=strcompress(string(science.vdisp_smm, format='(D10.1)'), /rem)+(science.vdisperr_smm le 0 ? '' : ' +/- '+strcompress(string(science.vdisperr_smm, format='(D10.1)'), /rem))+' km/s (SMM)'
        else: widget_control, widget_info(self.base, find_by_uname='curvdisp'), set_value=unknown
    endcase
end


pro sps_fit::getscience, files=files
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Initializing ...'

    common mask_in, mask_in
    case 1 of
        strmid(mask_in, 0, 4) eq '0024' or strmid(mask_in, 0, 3) eq 'rse': begin
            photfile = '/scr2/nichal/keck/deimos/Cl0024MS0451/Cl0024master.v7.fits.gz'
            vdispfile ='/scr2/nichal/keck/deimos/Cl0024MS0451/Cl0024_sigmas.txt'
            readcol, vdispfile, vd_objname, vd_vdisp, format='A,D', comment='#', /silent
        end
        strmid(mask_in, 0, 4) eq '0451' or strmid(mask_in,0,3) eq 'smm' or strmid(mask_in, 0, 4) eq 'mock': begin
            photfile = '/scr2/nichal/keck/deimos/Cl0024MS0451/MS0451master.v14.fits.gz'
            vdispfile = '/scr2/nichal/keck/deimos/Cl0024MS0451/MS0451_sigmas.txt'
            readcol, vdispfile, vd_objname, vd_vdisp, vd_vdisperr, format='A,D,D', comment='#', /silent
        end
        strmid(mask_in, 0, 6) eq 'Cl1604': begin
            photfile = getenv('UCI')+'Cl1604/catalogs/Cl1604.fits.gz'
        end
        else: message, mask_in+' is not associated with a cluster that I know.'
    endcase
    phot = mrdfits(photfile, 1, /silent)
    phot = phot[where(phot.dec gt -90 and phot.dec lt 90)]

    common npixcom, npix
    npix = 8192
    if strlowcase(mask_in) eq 'cl1604deimos' then npix = 10625
    if strlowcase(mask_in) eq 'cl1604lris' then npix = 2570

    observatory, 'keck', obs
    sciencefits = self.directory+(self.lowsn eq 1 ? 'sps_fit_lowsn.fits.gz' : 'sps_fit.fits.gz')
    if ~file_test(sciencefits) then begin
        if ~keyword_set(files) then message, 'You must specify the FILES keyword if a sps_fit.fits.gz file does not exist.'
        c = n_elements(files)
        masks = strarr(c)
        slits = strarr(c)
        objnames = strarr(c)
        for i=0,c-1 do begin
            basefile = file_basename(files[i])
            extensions = strsplit(basefile, '.', /extract)
            masks[i] = extensions[1]
            slits[i] = extensions[2]
            objnames[i] = extensions[3]
        endfor
        w = where(strmatch(objnames, '*serendip*') eq 0, cw)
        if cw gt 0 then begin
            masks = masks[w]
            slits = slits[w]
            objnames = objnames[w]
            files = files[w]
        endif else begin
            message, 'No spectra found.'
        endelse

        nspec = n_elements(objnames)
        self.nspec = nspec
        speclist = masks+' '+strtrim(string(slits), 2)+' '+objnames
        widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist

        scienceall = replicate({science}, nspec)
        wgood = bytarr(nspec)+1
        for i=0,nspec-1 do begin
            science = {science}
            science.skyfit = -1
            science.skylinemask = -1
            science.objname = objnames[i]
            science.mask = masks[i]
            science.slit = slits[i]
;;; Find failure files:
            failurefile = file_dirname(files[i], /mark_directory)+'failure.?.'+strtrim(string(science.slit, format='(I3)'), 2)+'.sav'
            failurefile = file_search(failurefile, count=nfail)
            failure1 = 0
            failure2 = 0
            for j=0,nfail-1 do begin
                basefailfile = file_basename(failurefile[j])
                extensions = strsplit(basefailfile, '.', /extract)
                chip = fix(extensions[2])
                case 1 of
                    chip le 4: failure1 = 1
                    chip ge 5: failure2 = 1
                endcase
            end            

            data1 = mrdfits(files[i], 1, hdr, /silent)
            if strlowcase(strmid(mask_in, 0, 6)) eq 'cl1604' then begin
                lambda = data1.lambda
                spec = data1.spec
                ivar = data1.ivar
                airmass = 1.0
            endif else begin
                airmass = sxpar(hdr, 'AIRMASS')
                ras = sxpar(hdr, 'RA_OBJ')
                decs = sxpar(hdr, 'DEC_OBJ')
                decs_mask = sxpar(hdr, 'DEC')
                decssplit = strsplit(decs, ':', /extract)
                if stregex(decssplit[0], '\*', /boolean) eq 1 then begin
                    decssplit_mask = strsplit(decs_mask, ':', /extract)
                    decssplit[0] = decssplit_mask[0]
                    decs = strjoin(decssplit, ':')
                endif
                if stregex(decssplit[1], '\*', /boolean) eq 1 then begin
                    decssplit_mask = strsplit(decs_mask, ':', /extract)
                    decssplit[0] = '-00'
                    decssplit[1] = '00'
                    decs = strjoin(decssplit, ':')
                endif
                get_coords, coords, instring=ras+'  '+decs
                science.ra = coords[0]*15.
                science.dec = coords[1]
                science.jdobs = double(sxpar(hdr, 'MJD-OBS')) + 2400000.5

                data2 = mrdfits(files[i], 2, /silent)            
                ;;; Flag regions of the spectrum with bad wavelength solution:
                if failure1 then data1.lambda -= 10000.0
                if failure2 then data2.lambda += 10000.0

                lambda = [data1.lambda, data2.lambda]
                spec = [data1.spec, data2.spec]
                ivar = [data1.ivar, data2.ivar]
                skyspec = [data1.skyspec, data2.skyspec] 
            endelse

            if strlowcase(strmid(mask_in, 0, 6)) eq 'cl1604' then begin
                match, strlowcase(strtrim(phot.mask, 2))+'_'+strtrim(phot.slit, 2), strlowcase(strtrim(science.mask, 2))+'_'+strtrim(science.slit, 2), w1, w2
            endif else begin
                spherematch, phot.ra, phot.dec, science.ra, science.dec, 1.0/3600., w1, w2
            endelse
            if w1[0] eq -1 then begin
            endif else begin
                science.z = phot[w1].z
                science.zquality = phot[w1].zquality
                if strlowcase(strmid(mask_in, 0, 6)) eq 'cl1604' then begin
                    science.objname = phot[w1].phot_id
                    objnames[i] = phot[w1].phot_id
                    science.zspec = 0.0
                    science.ra = phot[w1].ra
                    science.dec = phot[w1].dec
                    science.r = phot[w1].rmag
                    science.i = phot[w1].imag
                    science.f606w = phot[w1].f606wmag
                    science.f814w = phot[w1].f814wmag
                    science.logmstar = phot[w1].logMstar_SED_PFW
                    science.age = phot[w1].age_SED / 1d9
                endif else begin
                    science.zspec = science.z
                    science.zsource = phot[w1].zsource
                    science.b = phot[w1].b_auto
                    science.v = phot[w1].v_auto
                    science.r = phot[w1].r_auto
                    science.i = phot[w1].i_auto
                    science.j = phot[w1].j_auto
                    science.k = phot[w1].k_auto
                    science.f814w = phot[w1].f814w_auto
                    science.berr = phot[w1].b_auto_err
                    science.verr = phot[w1].v_auto_err
                    science.rerr = phot[w1].r_auto_err
                    science.ierr = phot[w1].i_auto_err
                    science.jerr = phot[w1].j_auto_err
                    science.kerr = phot[w1].k_auto_err
                    science.f814werr = phot[w1].f814w_auto_err
                endelse

                science.vdisp_smm = -999d
                science.vdisperr_smm = -999d
                if (size(vd_objname))[0] gt 0 then begin
                    match, strtrim(vd_objname, 2), strtrim(science.objname, 2), w1, w2, count=cmatch
                    if cmatch gt 0 then begin
                        science[w2].vdisp_smm = vd_vdisp[w1]
                        if strmid(mask_in, 0, 4) eq '0451' then begin
                            science[w2].vdisperr_smm = vd_vdisperr[w1]
                        endif else begin
                            science[w2].vdisperr_smm = -999d
                        endelse
                    endif
                endif
            endelse
            case 1 of
                science.b gt 10.0 and science.b lt 50.0 and science.v gt 10.0 and science.v lt 50.0: science.phot_color = 'BV'
                science.b gt 10.0 and science.b lt 50.0 and science.r gt 10.0 and science.r lt 50.0: science.phot_color = 'BR'
                science.v gt 10.0 and science.v lt 50.0 and science.i gt 10.0 and science.i lt 50.0: science.phot_color = 'VI'
                science.v gt 10.0 and science.v lt 50.0 and science.k gt 10.0 and science.k lt 50.0: science.phot_color = 'VK'
                science.v gt 10.0 and science.v lt 50.0 and science.j gt 10.0 and science.j lt 50.0: science.phot_color = 'VJ'
                science.j gt 10.0 and science.j lt 50.0 and science.k gt 10.0 and science.k lt 50.0: science.phot_color = 'JK'
                else: science.phot_color = 'BV'
            endcase

            science.spec1dfile = files[i]
            science.good = 1
            science.goodsky = 0
            w = where(science.age le 0, c)
            if c gt 0 then begin
                science[w].age = -999d
                science[w].ageerr = -999d
            endif
            science.feh = -999d
            science.feherr = -999d
            science.vdisp = -999d
            science.vdisperr = -999d
            science.alphafe = 0.0

            self.i = i

            w = where(ivar gt 0 and finite(ivar), civar)
            if civar gt 10 then if w[civar-1] ne n_elements(ivar)-1 then ivar[w[civar-11:civar-1]] = 0
            if self.lowsn eq 1 then begin
                common random, seed
                factor = 3.0
                spec[w] += factor*ivar[w]^(-0.5)*randomn(seed, civar, /double, /normal)
                ivar[w] /= (1.0 + factor^2.)
            endif
            n = n_elements(lambda)

            tell = mrdfits('/scr2/nichal/workspace2/telluric/deimos_telluric_1.0.fits', 1, /silent)
            wtell = n_elements(tell)-1
            tell = tell[wtell]
            ptr_free, self.tell
            self.tell = ptr_new(tell)

            t = (-1*ts_diff(lambda, 1))[0:n-2]
            wt = where(t le 0, ct)
            if ct gt 0 then begin
                message, 'Wavelength array for '+strtrim(objnames[i], 2)+' is not monotonic.  Omitting.', /info
                wgood[i] = 0
                continue
            endif

            science.lambda = lambda
            science.spec = spec
            science.ivar = ivar
            science.airmass = airmass
            if strlowcase(strmid(mask_in, 0, 6)) ne 'cl1604' then begin
                skyspec /= 40.0*median(skyspec)
                science.skyspec = skyspec
            endif

            ;self->skytweak, science
            self->specres, science
            self->telluric, science
            self->continuum, science
            if array_equal(science.continuum, replicate(-999, n_elements(science.continuum))) then begin
                message, 'Failed at finding continuum for '+strtrim(objnames[i], 2)+'.  Omitting.', /info
                wgood[i] = 0
                continue
            endif
            self->mask, science
            science.spscont = 1.0
            self->indices, science, /noredraw

            self->statusbox, science=science
            scienceall[i] = science
        endfor
        self.i = 0
        wgood = where(wgood eq 1, cgood)
        scienceall = scienceall[wgood]
        ptr_free, self.science
        self.science = ptr_new(scienceall)
        self.nspec = cgood
        self->writescience
        self->specres_mask, self.directory
        speclist = masks[wgood]+' '+strtrim(string(slits[wgood]), 2)+' '+objnames[wgood]
        widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist
        widget_control, widget_info(self.base, find_by_uname='mode'), set_value=4
        ;self->fit_all
        self->writescience
     endif else begin
        scienceall = mrdfits(sciencefits, 1, /silent)

        self.nspec = n_elements(scienceall)
        speclist = scienceall.mask+' '+strtrim(string(scienceall.slit), 2)+' '+scienceall.objname
        widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist

        tell = (mrdfits('/scr2/nichal/workspace2/telluric/deimos_telluric_1.0.fits', 1, /silent))
        wtell = n_elements(tell)-1
        tell = tell[wtell]
        ptr_free, self.tell
        self.tell = ptr_new(tell)
        ptr_free, self.science
        self.science = ptr_new(scienceall)
        widget_control, widget_info(self.base, find_by_uname='mode'), set_value=4
    endelse
end


pro sps_fit::writescience
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Writing to database ...'
    scienceall = *self.science
    sciencefits = self.directory+(self.lowsn eq 1 ? 'sps_fit_lowsn.fits' : 'sps_fit.fits')
    mwrfits, scienceall, sciencefits, /create, /silent
    spawn, 'gzip -f '+sciencefits
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end


pro sps_fit::initialize_directory, directory=directory
    common mask_in, mask_in
    newdirectory = '/scr2/nichal/workspace2/sps_fit/data/'+mask_in+'/'
    if ~file_test(newdirectory) then file_mkdir, newdirectory

    if strmid(directory, 0, 1, /reverse_offset) ne '/' then directory = directory + '/'
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Reading directory ...'

    countfiles:
    files = file_search(directory, 'spec1d*.{fits,fits.gz}', count=c)
    sciencefits = newdirectory+(self.lowsn eq 1 ? 'sps_fit_lowsn.fits.gz' : 'sps_fit.fits.gz')
    if c eq 0 then begin
        files = file_search(directory+'*/spec1d*.{fits,fits.gz}', count=c)
    endif
    if c eq 0 and ~file_test(sciencefits) then begin
        message, 'Unknown mask.'
    endif

    self.directory = newdirectory

    self->getscience, files=files
    self.i = 0
    science = (*self.science)[self.i]
    self->statusbox, science=science
    self->default_range
    self.lambdalim = [3300, 7000]
    self->newspec, increment=0
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end


; =============== INIT ================
function sps_fit::INIT, directory=directory, lowsn=lowsn
    common sps_spec, sps, spsz, spsage
    if (size(sps))[1] eq 0 then spsspec = sps_interp(0.0, 5.0)

    base = widget_base(/row, title='sps_fit', uvalue=self, mbar=menu, tab_mode=0, units=1)
    file_menu = widget_button(menu, value='File', /menu)
    wexit = widget_button(file_menu, value='Save', uvalue='save', uname='save')
    wexit = widget_button(file_menu, value='Exit', uvalue='exit', uname='exit')
    tools_menu = widget_button(menu, value='Tools', /menu)
    wdefaultrange = widget_button(tools_menu, value='Default Spectrum Settings', uname='default_range', uvalue='default_range')
    wdefault_cont = widget_button(tools_menu, value='Default Continuum Regions', uname='default_cont', uvalue='default_cont')
    wdefault_mask = widget_button(tools_menu, value='Default Pixel Mask', uname='default_mask', uvalue='default_mask')
    wdefault_goodspec = widget_button(tools_menu, value='Default Good Spectrum', uname='default_goodspec', uvalue='default_goodspec')
    wreprepare_all = widget_button(tools_menu, value='Reprepare All', uname='reprepare_all', uvalue='reprepare_all')
    wreprepare_all = widget_button(tools_menu, value='Fit All', uname='fit_all', uvalue='fit_all')

    wleft = widget_base(base, /column, uname='left')
    wright = widget_base(base, /column, uname='right')
    widget_control, /managed, base

    ; ------ LEFT -------
    wplotmode = widget_base(wleft, /column, /align_center, /frame)
    wplotradio = cw_bgroup(wplotmode, ['telluric cross-correlation', 'continuum fit', 'continuum division', 'sky line fit', 'rest frame'], /column, /exclusive, set_value=4, uname='mode', uvalue='mode', /no_release)
    wstep = widget_base(wleft, /row, /align_center)
    wbackward = widget_button(wstep, value='<---', uvalue='backward', uname='backward', tab_mode=1, xsize=75)
    wforward = widget_button(wstep, value='--->', uvalue='forward', uname='forward', tab_mode=1, xsize=75)
    wprepbase = widget_base(wleft, /row, /align_center)
    wreprepare = widget_button(wprepbase, value='Reprepare', uvalue='reprepare', uname='reprepare', tab_mode=1, xsize=75)
    wfit = widget_button(wprepbase, value='Fit', uvalue='fit', uname='fit', tab_mode=1, xsize=75)
    wfitmgfe = widget_button(wprepbase, value='Fit[Mg/Fe]', uvalue='fit[Mg/Fe]', uname='fit[Mg/Fe]', tab_mode=1, xsize=75)
    windicesbase = widget_base(wleft, /row, /align_center)
    windices = widget_button(windicesbase, value='Compute Indices', uvalue='indices', uname='indices', tab_mode=1, xsize=125)
    wgoodbase = widget_base(wleft, /column, /align_center)
    wgood = cw_bgroup(wgoodbase, ['good sky', 'good spectrum','good fit'], /nonexclusive, set_value=[0, 0, 0], uname='good', uvalue='good')

    wcurobj = widget_base(wleft, /column, /align_center, tab_mode=0, /frame)
    widbase = widget_base(wcurobj, /align_left, /row, xsize=235)
    widlabel = widget_label(widbase, value='object ID:', /align_right, uname='idlabel', xsize=65)
    wcurid = widget_label(widbase, value='     ', /align_left, uname='curid', uvalue='curid', xsize=180)
    wvbase = widget_base(wcurobj, /align_center, /row)
    wvlabel = widget_label(wvbase, value='B = ', /align_right, uname='maglabel', xsize=95)
    wcurmag = widget_label(wvbase, value='     ', /align_left, uname='curmag', uvalue='curmag', xsize=150)
    wbvbase = widget_base(wcurobj, /align_center, /row)
    wbvlabel = widget_label(wbvbase, value='B-V = ', /align_right, uname='collabel', xsize=95)
    wcurcol = widget_label(wbvbase, value='     ', /align_left, uname='curcol', uvalue='curcol', xsize=150)
    wagebase = widget_base(wcurobj, /align_center, /row)
    wagelabel = widget_label(wagebase, value='age = ', /align_right, uname='agelabel', xsize=95)
    wcurage = widget_label(wagebase, value='     ', /align_left, uname='curage', uvalue='curage', xsize=150)
    wmstarbase = widget_base(wcurobj, /align_center, /row)
    wmstarlabel = widget_label(wmstarbase, value='log M* = ', /align_right, uname='mstarlabel', xsize=95)
    wcurmstar = widget_label(wmstarbase, value='     ', /align_left, uname='curmstar', uvalue='curmstar', xsize=150)
    wfehbase = widget_base(wcurobj, /align_center, /row)
    wfehlabel = widget_label(wfehbase, value='[Fe/H] = ', /align_right, uname='fehlabel', xsize=95)
    wcurfeh = widget_label(wfehbase, value='     ', /align_left, uname='curfeh', uvalue='curfeh', xsize=150)
    wfebase = widget_base(wcurobj, /align_center, /row)
    wfelabel = widget_label(wfebase, value='Fe = ', /align_right, uname='felabel', xsize=95)
    wcurfe = widget_label(wfebase, value='     ', /align_left, uname='curfe', uvalue='curfe', xsize=150)
    wmgbase = widget_base(wcurobj, /align_center, /row)
    wmglabel = widget_label(wMgbase, value='Mg = ', /align_right, uname='mglabel', xsize=95)
    wcurMg = widget_label(wMgbase, value='     ', /align_left, uname='curmg', uvalue='curmg', xsize=150)
    wcahbase = widget_base(wcurobj, /align_center, /row)
    wcahlabel = widget_label(wcahbase, value='CaH = ', /align_right, uname='cahlabel', xsize=95)
    wcurcah = widget_label(wcahbase, value='     ', /align_left, uname='curcah', uvalue='curcah', xsize=150)
    wgbandbase = widget_base(wcurobj, /align_center, /row)
    wgbandlabel = widget_label(wgbandbase, value='Gband = ', /align_right, uname='gbandlabel', xsize=95)
    wcurgband = widget_label(wgbandbase, value='     ', /align_left, uname='curgband', uvalue='curgband', xsize=150)
    ;wafebase = widget_base(wcurobj, /align_center, /row)
    ;wafelabel = widget_label(wafebase, value='[a/Fe] = ', /align_right, uname='afelabel', xsize=95)
    ;wcurafe = widget_label(wafebase, value='     ', /align_left, uname='curafe', uvalue='curafe', xsize=150)
    wvdispbase = widget_base(wcurobj, /align_center, /row)
    wvdisplabel = widget_label(wvdispbase, value='sigma_v = ', /align_right, uname='vdisplabel', xsize=95)
    wcurvdisp = widget_label(wvdispbase, value='     ', /align_left, uname='curvdisp', uvalue='curvdisp', xsize=150)
    wsnbase = widget_base(wcurobj, /align_center, /row)
    wsnlabel = widget_label(wsnbase, value='SNR = ', /align_right, uname='snlabel', xsize=95)
    wcursn = widget_label(wsnbase, value='     ', /align_left, uname='cursn', uvalue='cursn', xsize=150)
    wzbase = widget_base(wcurobj, /align_center, /row)
    wzlabel = widget_label(wzbase, value='z = ', /align_right, uname='zlabel', xsize=95)
    wcurz = widget_label(wzbase, value='     ', /align_left, uname='curz', uvalue='curz', xsize=150)
    wzfitbase = widget_base(wcurobj, /align_center, /row)
    wzfitlabel = widget_label(wzfitbase, value='zfit = ', /align_right, uname='zfitlabel', xsize=95)
    wcurzfit = widget_label(wzfitbase, value='     ', /align_left, uname='curzfit', uvalue='curzfit', xsize=150)
    wzqualitybase = widget_base(wcurobj, /align_center, /row)
    wzqualitylabel = widget_label(wzqualitybase, value='zquality = ', /align_right, uname='zqualitylabel', xsize=95)
    wcurzquality = widget_label(wzqualitybase, value='     ', /align_left, uname='curzquality', uvalue='curzquality', xsize=150)
    wchisqbase = widget_base(wcurobj, /align_center, /row)
    wchisqlabel = widget_label(wchisqbase, value='chisq = ', /align_right, uname='chisqlabel', xsize=95)
    wcurchisq = widget_label(wchisqbase, value='     ', /align_left, uname='curchisq', uvalue='curchisq', xsize=150)
    wnloopbase = widget_base(wcurobj, /align_center, /row)
    wnlooplabel = widget_label(wnloopbase, value='nloop = ', /align_right, uname='nlooplabel', xsize=95)
    wcurnloop = widget_label(wnloopbase, value='     ', /align_left, uname='curnloop', uvalue='curnloop', xsize=150)

    wcursor = widget_base(wleft, /column, /align_right, tab_mode=0, /frame)
    wobswave = widget_label(wcursor, value='     ', /align_right, uname='obswave', uvalue='obswave', xsize=150)
    wrestwave = widget_label(wcursor, value='     ', /align_right, uname='restwave', uvalue='restwave', xsize=150)
    wyslit = widget_label(wcursor, value='     ', /align_right, uname='yslit', uvalue='yslit', xsize=150)

    ; ------ RIGHT -------
    wfile = widget_base(wright, /frame, /row, /align_left, tab_mode=1)
    wback = widget_button(wfile, value='Back', uvalue='back', uname='back', tab_mode=1)
    wfilelist = widget_combobox(wfile, uname='filelist', value='                 ', tab_mode=1, /dynamic_resize)
    wnext = widget_button(wfile, value='Next', uvalue='next', uname='next', tab_mode=1)
    wstatus = widget_text(wfile, xsize=108, value='Initializing ...', uname='status', uvalue='status', tab_mode=0)

    wspec = widget_base(wright, /frame, /column)
    wspecplot = widget_draw(wspec, xsize=1600, ysize=400, uname='spec', /button_events, keyboard_events=1)
    w2dplot = widget_draw(wspec, xsize=1600, ysize=300, uname='2d', /tracking_events, /motion_events)

    wspeccontrol = widget_base(wright, /row, /align_center, tab_mode=1)
    wycontrol = widget_base(wspeccontrol, /frame, /row)
    wylow = widget_text(wycontrol, xsize=8, /editable, uname='ylow', uvalue='ylow')
    wylabel = widget_label(wycontrol, value=' < y < ', /align_center, uname='ylabel')
    wyhigh = widget_text(wycontrol, xsize=8, /editable, uname='yhigh', uvalue='yhigh')
    wlambdacontrol = widget_base(wspeccontrol, /frame, /row)
    wblue = widget_button(wlambdacontrol, value='<-', uname='blue', uvalue='blue', /align_center)
    wlambdalow = widget_text(wlambdacontrol, xsize=8, /editable, uname='lambdalow', uvalue='lambdalow')
    wlambdalabel = widget_label(wlambdacontrol, value=' < l < ', /align_center, uname='lambdalabel', font='-urw-standard symbols l-medium-r-normal--0-0-0-0-p-0-adobe-symbol')
    wlambdahigh = widget_text(wlambdacontrol, xsize=8, /editable, uname='lambdahigh', uvalue='lambdahigh')
    wred = widget_button(wlambdacontrol, value='->', uname='red', uvalue='red', /align_center)

    widget_control, base, /realize
    self.base = base
    xmanager, 'sps_fit', self.base, /no_block, cleanup='sps_fit_cleanup'

    readcol,'/scr2/nichal/workspace2/telluric/telluric.mask', tellstart, tellend, format='D,D', /silent, comment='#'
    wbands = [1,2,3,4,5]
    tellstart = tellstart[wbands]
    tellend = tellend[wbands]
    ptr_free, self.linewaves, self.linewaves, self.linecolors, self.tellstart, self.tellend, self.tellthick
    self.linewaves = ptr_new([2798.0, 3646.00, 3727.425, 3750.15, 3770.63, 3797.90, 3835.39, 3868.71, 3888.65, 3889.05, 3933.663, 3967.41, 3968.468, 3970.07, 4101.76, 4305.05, 4340.47, 4861.33, 4958.92, 5006.84, 5167.321, 5172.684, 5183.604, 5875.67, 5889.951, 5895.924, 6300.30, 6548.03, 6562.80, 6583.41, 6678.152, 6716.47, 6730.85])
    self.linenames = ptr_new(['MgII', 'Hbreak', '[OII]', 'H12', 'H11', 'H10', 'H9', '[NeIII]', 'HeI', 'H8', 'CaH', '[NeIII]', 'CaK', 'He', 'Hd', 'CH', 'Hg', 'Hb', '[OIII]', '[OIII]', 'Mgb', 'Mgb', 'Mgb', 'HeI', 'NaD', 'NaD', '[OI]', '[NII]', 'Ha', '[NII]', 'HeI', '[SII]', '[SII]'])
    self.linecolors = ptr_new(['blue', 'black', 'blue', 'black', 'black', 'black', 'black', 'blue', 'blue', 'black', 'red', 'blue', 'red', 'black', 'black', 'red', 'black', 'black', 'blue', 'blue', 'red', 'red', 'red', 'blue', 'red', 'red', 'blue', 'blue', 'black', 'blue', 'blue', 'blue', 'blue'])
    self.tellstart = ptr_new(tellstart)
    self.tellend = ptr_new(tellend)
    self.tellthick = ptr_new([5, 2, 5, 2, 2])

    readcol, '/scr2/nichal/workspace2/sps_fit/lines.txt', linestart, lineend, linetype, format='D,D,A,X', /silent, comment='#'
    self.linestart = ptr_new(linestart)
    self.lineend = ptr_new(lineend)
    self.linetype = ptr_new(linetype)

    readcol,'/scr2/nichal/workspace2/sps_fit/lick_indices.txt',indnum,indbandstart,indbandend,bluecontstart,bluecontend,redcontstart,redcontend,junk,indname,format='I,D,D,D,D,D,D,I,A',comment='#',/silent
    self.indstart = ptr_new(indbandstart)
    self.indend   = ptr_new(indbandend)
;    self.indstart = ptr_new(indbandstart)
;    self.indend   = ptr_new(indbandend)
    self.indname  = ptr_new(indname)

    common random, seed
    seed = systime(1)
    
    self.i = -1
    self.lowsn = keyword_set(lowsn) ? 1 : 0
    self->initialize_directory, directory=directory
    return, 1
end

pro sps_fit__define 
    state = {sps_fit, $
             base:0L, $
             directory:'', $
             science:ptr_new(), $
             tell:ptr_new(), $
             lambdalim:[-100d, 100d], $
             ylim:[-100d, 100d], $
             divylim:[-100d, 100d], $
             skyylim:[-100d, 100d], $
             linewaves:ptr_new(), $
             linenames:ptr_new(), $
             linecolors:ptr_new(), $
             tellstart:ptr_new(), $
             tellend:ptr_new(), $
             tellthick:ptr_new(), $
             linestart:ptr_new(), $
             lineend:ptr_new(), $
             linetype:ptr_new(), $
             indstart:ptr_new(), $
             indend:ptr_new(), $
             indname:ptr_new(), $
             nspec:0L, $
             i:0L, $
             keystate:0, $
             lambda1:0d, $
             lowsn:0b}
end


pro science__define
    common npixcom, npix
    nsky = 200
    science = {science, $
               objname:'', $
               mask:'', $
               slit:0L, $
               ra:-999d, $
               dec:-999d, $
               jdobs:-999d, $
               lambda:dblarr(npix), $
               spec:dblarr(npix), $
               ivar:dblarr(npix), $
               skyspec:dblarr(npix), $
               telldiv:dblarr(npix), $
               telldivivar:dblarr(npix), $
               contmask:bytarr(npix), $
               continuum:dblarr(npix), $
               spscont:dblarr(npix), $
               contdiv:dblarr(npix), $
               contdivivar:dblarr(npix), $
               spsspec:dblarr(npix), $
               spsspecfull:dblarr(npix), $
               spscontfull:dblarr(npix), $
               spsspecmg:dblarr(npix), $
               spsspecfullmg:dblarr(npix), $
               spscontfullmg:dblarr(npix), $
               spsspecfe:dblarr(npix), $
               spsspecfullfe:dblarr(npix), $
               spscontfullfe:dblarr(npix), $
               fitmask:bytarr(npix), $
               fitMgmask:bytarr(npix), $
               fitFemask:bytarr(npix), $
               dlam:dblarr(npix), $
               resscale:-999d, $
               skyfit:[[dblarr(nsky)], [dblarr(nsky)], [dblarr(nsky)]], $
               skylinemask:lonarr(nsky), $
               goodsky:0, $
               airmass:-999d, $
               zobs:-999d, $
               zobserr:-999d, $
               zspec:-999d, $
               zfit:-999d, $
               z:-999d, $
               zquality:-999d, $
               zsource:0, $
               phot_color:'', $
               b:-999d, $
               v:-999d, $
               r:-999d, $
               i:-999d, $
               j:-999d, $
               k:-999d, $
               f606w:-999d, $
               f814w:-999d, $
               berr:-999d, $
               verr:-999d, $
               rerr:-999d, $
               ierr:-999d, $
               jerr:-999d, $
               kerr:-999d, $
               f606werr:-999d, $
               f814werr:-999d, $
               age:-999d, $
               ageerr:-999d, $
               mg:-999d, $
               mgerr:-999d, $
               fe:-999d, $
               feerr:-999d, $
               feh:-999d, $
               feherr:-999d, $
               alphafe:-999d, $
               alphafeerr:-999d, $
               logmstar:-999d, $
               b300:-999d, $
               vdisp:-999d, $
               vdisperr:-999d, $
               vdisp_smm:-999d, $
               vdisperr_smm:-999d, $
               vdisp_ppxf:-999d, $
               vdisperr_ppxf:-999d, $
               cah:-999d, $
               caherr:-999d, $
               gband:-999d, $
               gbanderr:-999d, $
               chisq:-999d, $
               chisqmg:-999d, $
               chisqfe:-999d, $
               sn:-999d, $
               spec1dfile:'', $
               good:0B,$
               goodfit:0B, $
               nloop:0}
end


pro sps_fit_ppxf, mask
    common mask_in, mask_in
    mask_in = mask
    
    directory = '/scr2/nichal/keck/deimos/Cl0024MS0451/enk_reduced/'+mask

    if strlowcase(mask) eq 'smm1a' then directory = '/scr2/nichal/keck/deimos/Cl0024MS0451/nicha_reduced/'+mask
    if strlowcase(mask) eq 'mock' then directory = '/scr2/nichal/workspace2/mock_data/mock'
    if strlowcase(mask) eq 'mock_es' then directory = '/scr2/nichal/workspace2/mock_data/mock_es'
    if strlowcase(mask) eq 'mockv3' then directory = '/scr2/nichal/workspace2/mock_data/mockv3'

    if strlowcase(mask) eq 'cl1604deimos' then directory = '/scr2/nichal/keck/deimos/Cl1604/deimos'
    if strlowcase(mask) eq 'cl1604lris' then directory = '/scr2/nichal/keck/deimos/Cl1604/lris'
    if ~file_test(directory) then message, 'Mask not found.'
    n = obj_new('sps_fit', directory=directory)
end
