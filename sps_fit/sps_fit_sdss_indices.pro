; =================
pro sps_iterproc, funcname, p, iter, fnorm, functargs=functargs, parinfo=pi, quiet=quiet, dof=dof
    common sps_iterproc, contiter
    if iter gt 1 then print, contiter, p[0], p[1], p[2],p[3], fnorm/dof, dof, format='(I3,2X,D6.3,1X,D5.2,2X,D6.1,2x,D6.3,1X,D8.3,1X,I4)'
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


pro sps_fit::fit, science, noredraw=noredraw, nostatusbar=nostatusbar
    common sps_spec, sps, spsz, spsage
    common sps_iterproc, contiter
    common get_sps, dlam, dataivar, datalam, wonfit
    makeplot=0

    if ~keyword_set(nostatusbar) then widget_control, widget_info(self.base, find_by_uname='status'), set_value='Fitting ...'
    restlambda = science.lambda / (1d + science.zspec)
    znow = science.zspec
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
    pi.value = double([0.0, 5.0, 200.0,znow])
    pi[0].limits = minmax(spsz)
    pi[1].limits = minmax(spsage)
    pi[2].limits = [0.0, 1000.0]
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
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, strtrim(science.objname, 2)+'  ('+strtrim(string(self.i+1, format='(I3)'), 2)+' / '+strtrim(string(self.nspec, format='(I3)'), 2)+')'
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, '  i Z/Z_sun   age sigma_v  redhift    chi^2  DOF'
    print, '--- ------- ----- ------- ---------  -------- ----'
    nloop = 0
    while abs(zdiff) gt 0.001 or abs(agediff) gt 0.001 or abs(vdispdiff) gt 0.001 or abs(redshfdiff) gt 0.0001 and nloop lt 30 do begin
        contiter++
        dlam = dlam_all
        dataivar = science.ivar*(median(science.spec))^2
        datalam  = science.lambda
        wonfit   = wontofit
        pars = mpfitfun('get_sps_obs_sdss', xmp, ymp, dymp, parinfo=pi, /nocatch, bestnorm=bestnorm, dof=dof, perror=perror, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=1000, status=status, yfit=ympfit, iterproc='sps_iterproc')
        
        zdiff = (pi[0].value-pars[0])/pi[0].value
        agediff = (pi[1].value-pars[1])/pi[1].value
        vdispdiff = (pi[2].value-pars[2])/pi[2].value
        redshfdiff = (pi[3].value-pars[3])/pi[3].value
        pi.value = pars
        nnan = where(finite(science.contdiv) and restlambda gt 3500 and restlambda lt 7400)
        wonfit  = nnan
        restlambda = science.lambda / (1d + pi[3].value)
        datalam = restlambda
        spsbestfitarr = get_sps_rest_sdss(restlambda[nnan], pi.value)
        spsbestfit = spsbestfitarr[*,0]
                                ;bkpt = slatec_splinefit(restlambda[nnan], science.contdiv[nnan] / spsbestfit, coeff, invvar=science.contdivivar[nnan], bkspace=165, upper=3, lower=3, /silent)
        
        bkpt = slatec_splinefit(restlambda[nnan], science.contdiv[nnan] / spsbestfit, coeff, bkspace=165, upper=3, lower=3, /silent)

        if bkpt[0] eq -1 then begin
            pi.value = [-999d, -999d, -999d,-999d]
            perror = [-999d, -999d, -999d, -999d]
            science.spsspec = -999d
            science.spscont = -999d ;nicha added this line
            break
        endif
        cont = slatec_bvalu(restlambda, bkpt, coeff)

        if makeplot eq 1 then makeloopplot,restlambda,science.contdiv,spsbestfit,cont,nnan,nloop,strtrim(string(self.i+1, format='(I3)'), 2),bestnorm/dof,pi.value,xmp/(1.+pi[3].value),ymp

        ymp = science.contdiv[won] / cont[won]
        dymp = (science.contdivivar[won])^(-0.5) / cont[won]
        science.spscont = cont ;nicha moved from below endwhile here
        nloop = nloop+1
        if nloop eq 30 then begin
           Print,'MAX OUT NLOOP!!!!'
        endif
     endwhile
 
    science.spsspec[nnan] = spsbestfit 
    science.spsspecfull[nnan] = spsbestfitarr[*,1]
    science.spscontfull[nnan] = spsbestfitarr[*,2]
    print, ' '
    
    done:
    science.nloop = nloop
    science.feh = pi[0].value
    science.feherr = perror[0]
    science.age = pi[1].value
    science.ageerr = perror[1]
    science.zfit = pi[3].value
    science.vdisp = pi[2].value
    science.vdisperr = perror[2]
    ;calculate chisq
    science.chisq = total((spsbestfit[won]-science.contdiv[won])^2*science.contdivivar[won])/float(n_elements(won))
    self->statusbox, science=science
    if ~keyword_set(noredraw) then begin
        self->redraw
    endif
end

pro sps_fit::fitmg, science, noredraw=noredraw, nostatusbar=nostatusbar
    common sps_spec, sps, spsz, spsage
    common sps_iterproc, contiter
    common get_sps, dlam, dataivar, datalam, wonfit

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
    zdiff = 1.0
    agediff = 1.0
    vdispdiff = 1.0
    redshfdiff = 1.0
    contiter = 0
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, strtrim(science.objname, 2)+'  ('+strtrim(string(self.i+1, format='(I3)'), 2)+' / '+strtrim(string(self.nspec, format='(I3)'), 2)+') Mg'
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, '  i Z/Z_sun   age sigma_v  redhift    chi^2  DOF'
    print, '--- ------- ----- ------- ---------  -------- ----'
    nloop = 0
    while abs(zdiff) gt 0.001 and nloop lt 1 do begin
        contiter++
        dlam = dlam_all
        dataivar = science.ivar*(median(science.spec))^2
        datalam  = science.lambda
        wonfit   = wontofit
        pars = mpfitfun('get_sps_obs_sdss', xmp, ymp, dymp, parinfo=pi, /nocatch, bestnorm=bestnorm, dof=dof, perror=perror, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=1000, status=status, yfit=ympfit, iterproc='sps_iterproc')
        
        zdiff = (pi[0].value-pars[0])/pi[0].value
        pi.value = pars
        nnan = where(finite(restlambda) and restlambda gt 3500 and restlambda lt 7400)
        wonfit  = nnan
        restlambda = science.lambda / (1d + pi[3].value)
        datalam = restlambda
        spsbestfitarr = get_sps_rest_sdss(restlambda[nnan], pi.value)
        spsbestfit = spsbestfitarr[*,0]
        bkpt = slatec_splinefit(restlambda[nnan], science.contdiv[nnan] / spsbestfit, coeff, invvar=science.contdivivar[nnan], bkspace=165, upper=3, lower=3, /silent)
        
        if bkpt[0] eq -1 then begin
            pi.value = [-999d, -999d, -999d,-999d]
            perror = [-999d, -999d, -999d, -999d]
            science.spsspec = -999d
            science.spscont = -999d ;nicha added this line
            break
        endif
        cont = slatec_bvalu(restlambda, bkpt, coeff)

        ymp = science.contdiv[won] / cont[won]
        dymp = (science.contdivivar[won])^(-0.5) / cont[won]
        science.spscont = cont ;nicha moved from below endwhile here
        nloop = nloop+1
        if nloop eq 30 then begin
           Print,'MAX OUT NLOOP!!!!'
        endif
     endwhile
 
    science.spsspecmg[nnan] = spsbestfit 
    science.spsspecfullmg[nnan] = spsbestfitarr[*,1]
    science.spscontfullmg[nnan] = spsbestfitarr[*,2]
    print, ' '
    
    done:
    science.mg = pi[0].value
    science.mgerr = perror[0]
 
    ;calculate chisq
    science.chisqmg = total((spsbestfit[won]-science.contdiv[won])^2*science.contdivivar[won])/float(n_elements(won))
    self->statusbox, science=science
    if ~keyword_set(noredraw) then begin
        self->redraw
    endif
end

pro sps_fit::fitfe, science, noredraw=noredraw, nostatusbar=nostatusbar
    common sps_spec, sps, spsz, spsage
    common sps_iterproc, contiter
    common get_sps, dlam, dataivar, datalam, wonfit

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
    while abs(zdiff) gt 0.001 and nloop lt 1 do begin
        contiter++
        dlam = dlam_all
        dataivar = science.ivar*(median(science.spec))^2
        datalam  = science.lambda
        wonfit   = wontofit
        pars = mpfitfun('get_sps_obs_sdss', xmp, ymp, dymp, parinfo=pi, /nocatch, bestnorm=bestnorm, dof=dof, perror=perror, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=1000, status=status, yfit=ympfit, iterproc='sps_iterproc')
        
        zdiff = (pi[0].value-pars[0])/pi[0].value
        pi.value = pars
        nnan = where(finite(restlambda) and restlambda gt 3500 and restlambda lt 7400)
        wonfit  = nnan
        restlambda = science.lambda / (1d + pi[3].value)
        datalam = restlambda
        spsbestfitarr = get_sps_rest_sdss(restlambda[nnan], pi.value)
        spsbestfit = spsbestfitarr[*,0]
        bkpt = slatec_splinefit(restlambda[nnan], science.contdiv[nnan] / spsbestfit, coeff, invvar=science.contdivivar[nnan], bkspace=165, upper=3, lower=3, /silent)
        
        if bkpt[0] eq -1 then begin
            pi.value = [-999d, -999d, -999d,-999d]
            perror = [-999d, -999d, -999d, -999d]
            science.spsspec = -999d
            science.spscont = -999d ;nicha added this line
            break
        endif
        cont = slatec_bvalu(restlambda, bkpt, coeff)

        ymp = science.contdiv[won] / cont[won]
        dymp = (science.contdivivar[won])^(-0.5) / cont[won]
        science.spscont = cont ;nicha moved from below endwhile here
        nloop = nloop+1
        if nloop eq 30 then begin
           Print,'MAX OUT NLOOP!!!!'
        endif
     endwhile
 
    science.spsspecfe[nnan] = spsbestfit 
    science.spsspecfullfe[nnan] = spsbestfitarr[*,1]
    science.spscontfullfe[nnan] = spsbestfitarr[*,2]
    print, ' '
    
    done:
    science.fe = pi[0].value
    science.feerr = perror[0]
 
    ;calculate chisq
    science.chisqmg = total((spsbestfit[won]-science.contdiv[won])^2*science.contdivivar[won])/float(n_elements(won))
    self->statusbox, science=science
    if ~keyword_set(noredraw) then begin
        self->redraw
    endif
 end

pro sps_fit::fit_all
    scienceall = *self.science
    curi = self.i
    for i=0,self.nspec-1 do begin
    ;for i=53,self.nspec-1 do begin
        self.i = i
        self->default_range
        science = scienceall[self.i]
        ;if science.good eq 0 then continue
        ;self->mask, science
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
    nmodes = 4
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode eq 3 and increment gt 0: mode = 0
        mode eq 0 and increment lt 0: mode = 3
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
    newi += increment
    if newi lt 0 then begin
       widget_control, widget_info(self.base, find_by_uname='status'), set_value='This is the first spectrum.'
       return
    endif
    if newi gt self.nspec-1 then begin
       widget_control, widget_info(self.base, find_by_uname='status'), set_value='This is the last spectrum.'
       return
    endif

    widget_control, widget_info(self.base, find_by_uname='filelist'), set_combobox_select=newi
    self.i = newi
    self->statusbox
    self.ylim = minmax((*self.science)[self.i].spec,/nan)
    self->skyylim
    ;self->default_range
    
    self.keystate = 0
    if ~keyword_set(noredraw) then self->redraw

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
        self.ylim = minmax((*self.science)[self.i].spec,/nan)
        self.ylim[1] *= 1.1
        self.divylim = [-1.0, 2.5]
        self->skyylim
        self.lambdalim = (minmax((*self.science)[self.i].lambda / (1d + (*self.science)[self.i].zspec),/nan) < 9100) > 2000
        self.lambdalim[0] >= 2000.
        self.lambdalim[1] <= 8938. / (1d + (*self.science)[self.i].zspec)
     endif
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode le 0: begin
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.ylim[0], format='(g8.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.ylim[1], format='(g8.2)'), /rem)
        end
        mode eq 2: begin
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.skyylim[0], format='(D5.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.skyylim[1], format='(D5.2)'), /rem)
        end
        else: begin
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.divylim[0], format='(D5.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.divylim[1], format='(D5.2)'), /rem)
        end
    endcase
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    zl = mode lt 3 ? (*self.science)[self.i].zspec : 0.0
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
    zl = mode lt 4 ? (*self.science)[self.i].zspec : 0.0
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
        mode le 0: self.ylim[0] = ylow
        mode eq 2: self.skyylim[0] = ylow
        else: self.divylim[0] = ylow
    endcase
    if ~keyword_set(noredraw) then self->redraw
end


pro sps_fit::yhigh, yhigh
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode le 0: self.ylim[1] = yhigh
        mode eq 2: self.skyylim[1] = yhigh
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



; ============== COMBOBOX  =============
pro sps_fit::handle_combobox, ev
    ;widget_control, ev.id, get_uvalue=uvalue
    ;widget_control, ev.top, get_uvalue=obj
    self.i = ev.index
    self->statusbox
    self.ylim = minmax((*self.science)[self.i].spec,/nan)
    self->skyylim
    self.keystate = 0
    self->redraw
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end


; ============= DRAW CLICK =============
pro sps_fit::handle_draw_click, ev
    click_coords = convert_coord(ev.x, ev.y, /device, /to_data)
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    if mode lt 3 then click_coords /= 1d + (*self.science)[self.i].zspec
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
            zl = mode lt 3 ? (*self.science)[self.i].zspec : 0.0
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
            zl = mode lt 3 ? (*self.science)[self.i].zspec : 0.0
            widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
            self->lambdalow, llownew, /noredraw
            self->lambdahigh, lhighnew
        end
        2: begin
            widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
            case mode of
                mode le 0: ylim = self.ylim
                mode eq 2: ylim = self.skyylim
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
                mode le 0: begin
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

        0: begin        ;continuum fit
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

        1: begin        ;continuum division
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

        2: begin        ;sky line fit
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

        3: begin        ;pixel mask
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
                            w = where(scienceall[self.i].lambda/(1d + scienceall[self.i].zspec) gt lambda1 and scienceall[self.i].lambda/(1d + scienceall[self.i].zspec) lt lambda2, c)
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
                            w = where(scienceall[self.i].lambda/(1d + scienceall[self.i].zspec) gt lambda1 and scienceall[self.i].lambda/(1d + scienceall[self.i].zspec) lt lambda2, c)
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
       z  = science.zspec
       if sn gt 6. and z gt 0.0 then science.good = 1 else science.good = 0
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
    update_else = 1

    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Repreparing ...'
    scienceall = *self.science
    science = scienceall[self.i]
    contmask = science.contmask

    n = n_elements(science.lambda)
    wcont = where(contmask[3:n-4] eq 1)+3
    wcont = wcont[where(finite(science.spec[wcont]) and finite(science.continuum[wcont]) and science.continuum[wcont] ne 0)]
    dev = abs((science.spec[wcont] - science.continuum[wcont]) / science.continuum[wcont])
    avgdev = mean(dev)
    w = where(dev lt 3.0*avgdev, c)
    if c gt 0 then science.sn = 1.0/mean(dev[w])

    if update_else eq 1 then begin
        ;self->skytweak, science
        ;science.skylinemask = -1
        ;science.contmask = 0
        self->specres, science, /goodoverride
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
    w = where(~finite(skyspec), c)
    if c gt 0 then skyspec[w] = 0
    skyspec = skyspec/max(skyspec)
    medskyspec = median(skyspec)

    deriv1skyspec = deriv(lambda, skyspec)
    deriv2skyspec = deriv(lambda, deriv1skyspec)

    thresh = 1.
    nlines = 1000
    while nlines gt 200 do begin
        w = where(abs(deriv1skyspec) lt 0.2 and deriv2skyspec lt -0.005 and skyspec gt thresh*medskyspec)
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
        contmask[n-3:n-1] = 0
    endif

    satbands = [6870, 7650]
    niter = 5

    wwhole = lindgen(n)
    spec = science.spec
    contmask(where(finite(spec) eq 0)) = 0
    contmask[wwhole[0:3]] = 0
    won = where(contmask[wwhole] eq 1, complement=woff, con) + wwhole[0]
    woff += wwhole[0]
    case 1 of
       spline: begin
          bkpt = slatec_splinefit(lambda[won], spec[won], coeff, invvar=science.ivar[won], bkspace=165, upper=5, lower=1.5,mask=mask, /silent,/everyn)
          if bkpt[0] eq -1 then return
          cont = slatec_bvalu(lambda[wwhole], bkpt, coeff)
       end
       poly: begin
          degree = 12
          norm = median(spec[won])
          a = [norm, replicate(0.0, degree-1)]
          p = lmfit(lambda[won], spec[won], a, measure_errors=(science.ivar[won])^(-0.5), /double, function_name='legendre_poly')
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
             cont = smooth_gauss_wrapper(lambda[wcont], spec[wcont], lambda[wwhole], 30.0, ivar1=science.ivar[wcont])
             wcont = where(abs(spec[ww]-cont[ww-wwhole[0]]) lt (science.ivar[ww])^(-0.5) and contmask[ww] eq 1, ccont) + ww[0]
             
             if array_equal(wcont, wcontold) then break
          endfor
       end
    endcase
    
    science.contmask = contmask
    science.continuum = cont
    science.contdiv = science.spec/cont
    science.contdivivar = science.ivar*cont^2.

    wcont = where(contmask[3:n-4] eq 1)+3
    wcont = wcont[where(finite(science.spec[wcont]) and finite(science.continuum[wcont]) and science.continuum[wcont] ne 0)]
    dev = abs((science.spec[wcont] - science.continuum[wcont]) / science.continuum[wcont])
    avgdev = mean(dev)
    w = where(dev lt 3.0*avgdev, c)
    if c gt 0 then science.sn = 1.0/mean(dev[w])
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
        mask = bytarr(n_elements(science.lambda))+1
        wem = where(linetype eq 'e', cemlines) ;where emission lines are
        for i=0,cemlines-1 do begin
           w=where(science.lambda/(1.+science.zspec) gt linestart(wem[i]) and science.lambda/(1.+science.zspec) lt lineend(wem[i]),cw)
           if cw gt 0 then mask[w]=0
        endfor
        w = where(science.ivar le 0 or ~finite(science.ivar) or science.lambda/(1d + science.zspec) lt 3650 or science.lambda/(1.+science.zspec) gt 7400, cw)
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

        0: begin                ;continuum fit
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

            plot, science.lambda, science.spec, xrange=self.lambdalim * (1d + science.zspec), yrange=self.ylim, xstyle=5, ystyle=5, background=fsc_color('white'), color=fsc_color('black'), /nodata
            for i=0,cstart-1 do begin
                x = ([science.lambda[wstart[i]], science.lambda[wstart[i]], science.lambda[wend[i]], science.lambda[wend[i]]] > (self.lambdalim[0] * (1d + science.zspec))) < (self.lambdalim[1] * (1d + science.zspec))
                y = [self.ylim[0], self.ylim[1], self.ylim[1], self.ylim[0]]
                polyfill, x, y, color=fsc_color('light cyan')
            endfor
            oplot, science.lambda, science.spec, color=fsc_color('black')
            oplot, science.lambda, science.continuum, color=fsc_color('green')
            if science.feh lt 3. then begin
               normfactor = median(science.spec)/median(science.spsspecfull)
               oplot, science.lambda,science.spsspecfull*normfactor,color=fsc_color('red')            
               oplot, science.lambda,science.spscontfull*normfactor,color=fsc_color('orange')            
            endif
            if science.fe lt 3. then begin
               normfactor = median(science.spec)/median(science.spsspecfullfe)
               oplot, science.lambda,science.spsspecfullfe*normfactor,color=fsc_color('darkgreen')            
               oplot, science.lambda,science.spscontfullfe*normfactor,color=fsc_color('darkgreen')            
            endif
            if science.mg lt 3. then begin
               normfactor = median(science.spec)/median(science.spsspecfullmg)
               oplot, science.lambda,science.spsspecfullmg*normfactor,color=fsc_color('blue')            
               oplot, science.lambda,science.spscontfullmg*normfactor,color=fsc_color('blue')            
            endif
            
            plot, science.lambda, science.spec, xrange=self.lambdalim * (1d + science.zspec), yrange=self.ylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6flux (e!E-!N/hr)!3', /nodata, /noerase
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

        1: begin                ;continuum division
            plot, science.lambda, science.contdiv, xrange=self.lambdalim * (1d + science.zspec), yrange=self.divylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6flux (normalized)!3', /nodata
            oplot, self.lambdalim, [1.0, 1.0], color=fsc_color('orange')
            oplot, science.lambda, science.contdiv, color=fsc_color('black')
            n = n_elements(*self.tellstart)
            for i=0,n-1 do begin
                oplot, [(*self.tellstart)[i], (*self.tellend)[i]], 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(*self.tellthick)[i]
            endfor
        end

        2: begin        ;sky line fit
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

        3: begin        ;rest frame
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
    zl = mode lt 3 ? (*self.science)[self.i].zspec : 0.0
    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)

    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end


pro sps_fit::statusbox, science=science
    if ~keyword_set(science) then science = (*self.science)[self.i]
    
    unknown = '???'
    widget_control, widget_info(self.base, find_by_uname='good'), set_value=[science.goodsky, science.good, science.goodfit]
    widget_control, widget_info(self.base, find_by_uname='curid'), set_value=strtrim(string(fix(science.plate)), 2)+'/'+strtrim(string(fix(science.mjd)), 2)+'/'+strtrim(string(fix(science.fiber)), 2)+'/'+' ('+strcompress(self.i+1, /rem)+' / '+strcompress(self.nspec, /rem)+')'
    widget_control, widget_info(self.base, find_by_uname='curz'), set_value=strcompress(string(science.zspec, format='(D5.3)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='curzfit'), set_value=strcompress(string(science.zfit, format='(D5.3)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='curchisq'), set_value=strcompress(string(science.chisq, format='(D4.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='curnloop'), set_value=strcompress(string(science.nloop, format='(D4.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='cursn'), set_value=science.sn gt 0 ? strcompress(string(science.sn, format='(D10.1)'), /rem) : unknown
    widget_control, widget_info(self.base, find_by_uname='curage'), set_value=science.age gt -100 ? strcompress(string(science.age, format='(D10.2)'), /rem)+(science.ageerr le 0 ? '' : ' +/- '+strcompress(string(science.ageerr, format='(D10.2)'), /rem))+' Gyr' : unknown
    widget_control, widget_info(self.base, find_by_uname='curmstar'), set_value=science.logmstar gt 0 ? strcompress(string(science.logmstar, format='(D10.2)'), /rem) : unknown
    widget_control, widget_info(self.base, find_by_uname='curfeh'), set_value=science.feh gt -100 ? strcompress(string(science.feh, format='(D10.2)'), /rem)+(science.feherr le 0 ? '' : ' +/- '+strcompress(string(science.feherr, format='(D10.2)'), /rem)) : unknown
    widget_control, widget_info(self.base, find_by_uname='curfe'), set_value=science.fe gt -100 ? strcompress(string(science.fe, format='(D10.2)'), /rem)+(science.feerr le 0 ? '' : ' +/- '+strcompress(string(science.feerr, format='(D10.2)'), /rem)) : unknown
    widget_control, widget_info(self.base, find_by_uname='curmg'), set_value=science.mg gt -100 ? strcompress(string(science.mg, format='(D10.2)'), /rem)+(science.mgerr le 0 ? '' : ' +/- '+strcompress(string(science.mgerr, format='(D10.2)'), /rem)) : unknown
    widget_control, widget_info(self.base, find_by_uname='curcah'), set_value=science.caherr gt 0 ? strcompress(string(science.cah, format='(D10.2)'), /rem)+(science.caherr le 0 ? '' : ' +/- '+strcompress(string(science.caherr, format='(D10.2)'), /rem)) : unknown
    widget_control, widget_info(self.base, find_by_uname='curvdisp'), set_value=strcompress(string(science.vdisp, format='(D10.1)'), /rem)+(science.vdisperr le 0 ? '' : ' +/- '+strcompress(string(science.vdisperr, format='(D10.1)'), /rem))+' km/s'

end


pro sps_fit::getscience, files=files
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Initializing ...'
    common mask_in, mask_in
    common npixcom, npix
    npix = 3700
    observatory, 'apo', obs
    sciencefits = self.directory+(self.lowsn eq 1 ? 'sps_fit_lowsn.fits.gz' : 'sps_fit.fits.gz')
    if ~file_test(sciencefits) then begin
        if ~keyword_set(files) then message, 'You must specify the FILES keyword if a sps_fit.fits.gz file does not exist.'
        c = n_elements(files)
        plates = strarr(c)
        MJDs = strarr(c)
        fibers = strarr(c)
        for i=0,c-1 do begin
            basefile = file_basename(files[i])
            extensions = strsplit(basefile, '-.', /extract)
            plates[i] = extensions[1]
            MJDs[i] = extensions[2]
            fibers[i] = extensions[3]
        endfor
        nspec = n_elements(fibers)
        self.nspec = nspec
        speclist = plates+' '+MJDs+' '+fibers
        widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist

        scienceall = replicate({science}, nspec)
        wgood = bytarr(nspec)+1
        for i=0,nspec-1 do begin
            science = {science}
            science.skyfit = -1
            science.skylinemask = -1
            science.mask = mask_in
            science.fiber = fibers[i]
            science.plate = plates[i]
            science.MJD = MJDs[i]

            data1 = mrdfits(files[i], 1,/silent)
            hdr   = mrdfits(files[i], 2,/silent)
               
            science.objname = hdr.specobjid
            science.zspec = hdr.z
           ;fix number of elements to match npix
            nlamb = n_elements(data1.loglam)
            if nlamb gt npix then data1=data1[0:npix-1]
            if nlamb lt npix then begin
               sdss_nan = data1[0]
               for ntag=0,n_elements(tag_names(data1))-1 do sdss_nan.(ntag) = 1./0.
               nanarr   = replicate(sdss_nan,npix-nlamb)
               nanarr.ivar = 0
               data1 = [data1,nanarr]
            endif
            lambda    = 10.^data1.loglam
            spec      = data1.flux
            ivar      = data1.ivar
            
            science.spec1dfile = files[i]
            science.good = 1
            science.goodsky = 1
            goodspec_tag = where(data1.AND_Mask eq 0,cgoodspec)
            if cgoodspec ne 0 then science.fitmask(goodspec_tag) = 1
            science.fitmask(where(~finite(science.spec))) = 0

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
      
            n = n_elements(lambda)

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
            science.skyspec = data1.sky
            science.sdssmodel = data1.model

            self->specres, science
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
        objlist = scienceall.objname
        speclist = mask_in+' '+strtrim(string(objlist[wgood]), 2)
        widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist
        widget_control, widget_info(self.base, find_by_uname='mode'), set_value=3
        self->writescience
     endif else begin
        scienceall = mrdfits(sciencefits, 1, /silent)
        self.nspec = n_elements(scienceall)
        speclist = mask_in+' '+strtrim(string(scienceall.objname), 2)
        widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist

        ptr_free, self.science
        self.science = ptr_new(scienceall)
        widget_control, widget_info(self.base, find_by_uname='mode'), set_value=3
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
    newdirectory = '/scr2/nichal/workspace2/sps_fit/data/indices/'+mask_in+'/'
    if ~file_test(newdirectory) then file_mkdir, newdirectory

    if strmid(directory, 0, 1, /reverse_offset) ne '/' then directory = directory + '/'
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Reading directory ...'

    countfiles:
    files = file_search(directory, 'spec-*.{fits,fits.gz}', count=c)
    sciencefits = newdirectory+(self.lowsn eq 1 ? 'sps_fit_lowsn.fits.gz' : 'sps_fit.fits.gz')
    if c eq 0 then begin
        files = file_search(directory+'*/spec-*.{fits,fits.gz}', count=c)
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

    base = widget_base(/row, title='sps_fit', uvalue=self, mbar=menu, tab_mode=0, units=1);,xsize=1000,ysize=600)
    ;------Top menu--------
    file_menu = widget_button(menu, value='File', /menu)
    wexit = widget_button(file_menu, value='Save', uvalue='save', uname='save')
    wexit = widget_button(file_menu, value='Exit', uvalue='exit', uname='exit')
    tools_menu = widget_button(menu, value='Tools', /menu)
    wdefaultrange = widget_button(tools_menu, value='Default Spectrum Settings', uname='default_range', uvalue='default_range')
    wdefault_cont = widget_button(tools_menu, value='Default Continuum Regions', uname='default_cont', uvalue='default_cont')
    wdefault_mask = widget_button(tools_menu, value='Default Pixel Mask', uname='default_mask', uvalue='default_mask')
    wdefault_maskall = widget_button(tools_menu, value='Default Pixel Mask All', uname='default_maskall', uvalue='default_maskall')
    wdefault_goodspec = widget_button(tools_menu, value='Default Good Spectrum', uname='default_goodspec', uvalue='default_goodspec')
    wreprepare_all = widget_button(tools_menu, value='Reprepare All', uname='reprepare_all', uvalue='reprepare_all')
    wreprepare_all = widget_button(tools_menu, value='Fit All', uname='fit_all', uvalue='fit_all')

    wleft = widget_base(base, /column, uname='left')
    wright = widget_base(base, /column, uname='right')
    widget_control, /managed, base

    ; ------ LEFT -------
    wplotmode = widget_base(wleft, /column, /align_center, /frame)
    wplotradio = cw_bgroup(wplotmode, ['continuum fit', 'continuum division', 'sky line fit', 'rest frame'], /column, /exclusive, set_value=3, uname='mode', uvalue='mode', /no_release)
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
    wchisqbase = widget_base(wcurobj, /align_center, /row)
    wchisqlabel = widget_label(wchisqbase, value='chisq = ', /align_right, uname='chisqlabel', xsize=95)
    wcurchisq = widget_label(wchisqbase, value='     ', /align_left, uname='curchisq', uvalue='curchisq', xsize=150)
    wnloopbase = widget_base(wcurobj, /align_center, /row)
    wnlooplabel = widget_label(wnloopbase, value='nloop = ', /align_right, uname='nlooplabel', xsize=95)
    wcurnloop = widget_label(wnloopbase, value='     ', /align_left, uname='curnloop', uvalue='curnloop', xsize=150)

    ; ------ RIGHT -------
    wfile = widget_base(wright, /frame, /row, /align_left, tab_mode=1)
    wback = widget_button(wfile, value='Back', uvalue='back', uname='back', tab_mode=1)
    wfilelist = widget_combobox(wfile, uname='filelist', value='                 ', tab_mode=1, /dynamic_resize)
    wnext = widget_button(wfile, value='Next', uvalue='next', uname='next', tab_mode=1)
    wstatus = widget_text(wfile, xsize=108, value='Initializing ...', uname='status', uvalue='status', tab_mode=0)

    wspec = widget_base(wright, /frame, /column)
    wspecplot = widget_draw(wspec, xsize=1600, ysize=700, uname='spec', /button_events, keyboard_events=1)

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
    ptr_free, self.linewaves, self.linewaves, self.linecolors, self.tellstart, self.tellend, self.tellthick, self.indstart, self.indend, self.indname
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
               MJD:-999d, $
               plate:-999d, $
               fiber:-999d, $
               lambda:dblarr(npix), $
               spec:dblarr(npix), $
               ivar:dblarr(npix), $
               skyspec:dblarr(npix), $
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
               sdssmodel:dblarr(npix), $
               fitmask:bytarr(npix), $
               fitMgmask:bytarr(npix), $
               fitFemask:bytarr(npix), $
               dlam:dblarr(npix), $
               skyfit:[[dblarr(nsky)], [dblarr(nsky)], [dblarr(nsky)]], $
               skylinemask:lonarr(nsky), $
               goodsky:1B, $
               nloop:-999d, $
               zspec:-999d, $
               zfit:-999d, $
               zsource:0, $
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
               vdisp:-999d, $
               vdisperr:-999d, $
               cah:-999d, $
               caherr:-999d, $
               gband:-999d, $
               gbanderr:-999d, $
               chisq:-999d, $
               chisqmg:-999d, $
               chisqfe:-999d, $
               sn:-999d, $
               spec1dfile:'', $
               good:1B,$
               goodfit:1B}
end


pro sps_fit_sdss_indices, mask
    common mask_in, mask_in
    mask_in = mask

    directory = '/scr2/nichal/workspace2/SDSSdata/'+mask
    if ~file_test(directory) then message, 'Mask not found.'
    n = obj_new('sps_fit', directory=directory)
end
