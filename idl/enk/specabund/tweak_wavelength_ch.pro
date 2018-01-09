function tweak_wavelength_ch, science, mds=mds
    if science.teff lt 3500 or science.teff gt 6400 or science.logg lt 0.0 or science.logg gt 4.0 or science.feh lt -4.0 or science.feh gt 0.0 or science.alphafe lt -0.8 or science.alphafe gt 1.2 or science.cfe lt -2.4 or science.cfe gt 3.5 then begin
        message, 'Atmospheric paramters are outside of the grid boundary.', /info
        return, replicate(0d, n_elements(science.lambda))
    endif

    common mooglambdacom_ch, mooglambda_ch
    common hashtable_ch, ht_ch, nmoog_ch
    moogspec_ch = interp_moog_ch(science.teff, science.logg, science.feh, science.alphafe, science.cfe, mds=mds)

    degree = 2
    errflag = 0

    lambda = science.lambda
    contdiv = science.contdiv
    contdivivar = science.contdivivar

    tmp = minmax(lambda)

    temp_flux = moogspec_ch
    temp_wave = mooglambda_ch * (1d + science.zrest)

;----restores template sky spectrum from HiRes: 
;----wavelength array in temp_wave, flux in temp_flux

;----smooth the HiRes spectrum with a gaussian to make similar to
;----DEIMOS resolution

    sigma = 9
    halfwidth = 15
    kernel = findgen(2*halfwidth+1) - halfwidth
    kernel = exp(-kernel^2 / 2 / sigma^2)
    kernel = kernel / total(kernel)

    temp_flux = convol(temp_flux, kernel, /center, /edge_wrap)

    sizex = n_elements(lambda)

    whgood = where(finite(contdivivar) and contdivivar gt 0, goodct)


;---set up 0.2A wavelength grid
    if goodct gt 2 then begin 
        minlambda = min(lambda[whgood], max=max1) > min(temp_wave, max=max2)
        maxlambda = max1 < max2
    endif else begin
        minlambda=0
        maxlambda=0
    endelse

    dlam = 0.2  ; size of shifts in cross correlation
    oversample = 0.2/dlam

    nlag = 4*oversample*2 + 1   ; number of shifts in cross correlation

    minlambda = (minlambda+1) < 1.5E4
    maxlambda = maxlambda-1 > 500      ; kill off spurious end effects


    npts = floor((maxlambda-minlambda)/dlam)
    fullwave = findgen(npts>1)*dlam + minlambda

;---interpolate observed spectrum onto grid
    iwave = where((lambda ge minlambda) and (lambda le maxlambda), ct)
    if ct gt 0 then fullflux = interpol(contdiv[iwave], lambda[iwave], fullwave) else fullflux = 0.

;---test whether bspline is useful:
    if total(fullflux) gt 0. then begin

;---interpolate template onto wavelength grid
        iwave = where((temp_wave ge minlambda) and (temp_wave le maxlambda))
        stemp_flux = interpol(temp_flux[iwave], temp_wave[iwave], fullwave)


        window = 100.                ; 100A windows
        nregions = fix(2*(maxlambda-minlambda) / window)
        
        color = maxlambda gt 8500. ? 'R' : 'B'

        ;if color eq 'B' then begin 
        ;    sig_thresh = 400.         ; minimum signal required in stemp_flux over window
        ;    peakthresh = 1.5          ; minimum height for a significant peak
        ;endif
        ;if color eq 'R' then begin 
        ;    sig_thresh = 400.           ; minimum signal required in stemp_wave over window
        ;    peakthresh = 10.           ;minimum signal required for a significant peak
        ;endif
        lag = indgen(nlag) - nlag/2  ;shift values for x-correlation
        wavlag = lag*dlam
        npoly = 3
        cc = fltarr(nlag, nregions)     ; cross-correlation values for lag
        cc_err = cc                     ;  estimated error in cross-correlation

        shiftfit = fltarr(npoly, nregions) ;polynomial coefficients for fit to cc
        fiterr = shiftfit                  ; errors on coefficients

        shift = fltarr(nregions)     ;peak of fit
        shifterr = shift             ;error in peak finding

        dofit = shift                 ;which regions' shifts have enough signal to trust
        fitlam = shift                ;central wavelength of fit regions
        fitpix = shift                ; central pixel of fit regions
        tmed = shift                  ;median of each window

        for i=0,nregions-1 do begin
;---march across spectrum in 100A windows, overlapping by 50A
            minwave = minlambda + i*window/2. + nlag/2*dlam + dlam
            maxwave = minwave + window
            if maxwave ge maxlambda then begin
                minwave = maxlambda - window - nlag/2*dlam - dlam
                maxwave = maxlambda
            endif

            iwin = where((fullwave ge minwave) and (fullwave le maxwave), cwin)
            if cwin lt 50 then begin
                shifterr[i] =  1d10
                dofit[i] = 0
            endif else begin
                fitpix[i] = min(iwin) + (max(iwin) - min(iwin))/2.
                wave = fullwave[iwin]
                fitlam[i] = fullwave[fitpix[i]]
                flux = fullflux[iwin]
                tflux = stemp_flux[iwin]
                djs_iterstat, tflux, median=tmp
                tmed[i] = tmp
                
;----does the input bspline exist here?
                if total(flux) ne 0. then begin    
;----is there enough signal to do cross-correlation?
                    ;if total(tflux-tmed[i]) ge sig_thresh then begin
                    

;----do cross-correlation, fit with a polynomial, find max
                    cc[*,i] = c_correlate(flux, tflux, lag)
                    
                    ;wpeak = where(tflux gt peakthresh) ;where there's a significant peak
                    ;if wpeak[0] ne -1 then error = .2/(sqrt(total(tflux[wpeak])/2.)) else error = 1.

                    ;if wpeak[0] ne -1 then print,total(tflux-tmed[i]),((total(tflux[wpeak])/2.))
                    ;cc_err[*,i] = replicate(error, nlag)

; ----only fit near the xcorr peak
                    maxcorr = max(cc[*,i], maxpix)
                    ; number of pix away from peak to go
                    nout = 2
                    mintofit = (maxpix-nout*oversample) > 0
                    if maxpix lt n_elements(wavlag) - nout*oversample - 1 then maxtofit = mintofit + 2*nout*oversample else begin
                        maxtofit = (maxpix + nout*oversample) < (n_elements(wavlag)-1)
                        mintofit = maxtofit - 2*nout*oversample
                    endelse

                    ;print, minwave,  maxwave,  cc_err[0]
                    ;shiftfit[*, i] =  poly_fit(wavlag[mintofit:maxtofit], cc[mintofit:maxtofit, i], npoly-1, measure_errors=cc_err[mintofit:maxtofit,i], sigma=err)
                    shiftfit[*, i] =  poly_fit(wavlag[mintofit:maxtofit], cc[mintofit:maxtofit, i], npoly-1, sigma=err)
                    fiterr[*,i] = err
                    shift[i] = -0.5*shiftfit[1, i]/shiftfit[2, i]
                    
                    shifterr[i] = (abs(shift[i])>0.005)*sqrt((fiterr[1,i]/shiftfit[1,i])^2 + (fiterr[2,i]/shiftfit[2,i])^2)
;---deweight if we're in the noisy part of the template.
                    ;if minwave ge 8950. then shifterr[i] = sqrt(5.)*shifterr[i]        
                    
                    dofit[i] = 1 
                    ;endif else begin
                    ;    shifterr[i] =  1d10
                    ;    dofit[i] = 0
                    ;endelse
                endif else begin
                    shifterr[i] =  1d10
                    dofit[i] = 0
                endelse
            endelse
        endfor

        whfit =  where(dofit)
        nfit = n_elements(whfit)
;whfit =  whfit[1:nfit-2];drop first and last points.
        fitlam = fitlam[whfit]
        shift = shift[whfit]
        shifterr = shifterr[whfit]
        fitpix = fitpix[whfit]
;  cent_row = floor(sizey/2.)
;  cent_wave = wave2d[*, cent_row]
        pix = findgen(sizex)

; find interpolated pixel numbers corresponding 
; to fit wavelengths
        fitpix2 = fitlam
        fitpix2 = interpol(pix, lambda, fitlam) 

; the new wavelengths for the fit pixels:
        shift_wave = shift ;fitlam+shift

        npix = n_elements(fullwave)
        npix2 = sizex

        shifterr=sqrt(shifterr^2+(median(shifterr) > 0.0025)^2)

; do new fits for regular grid and input wavelength array
        xx = fitpix2/(npix2/2.) -1

        pxx2 =  fitpix2

; fit vs. wavelength
;for degree = 2,3 do begin
        wave_fit2 = svdfit(xx, shift_wave, degree, /double, /legendre, yfit=pxx2, measure_errors = shifterr, chisq=chisq,sigma=sigma)
;  print,degree,chisq,djsig(shift_wave-pxx2)
;endfor

        chisq2 = chisq
        pxx =  fitpix



; fit vs. pixel in wavelength array
; this fit generally looks better, but makes less sense for global fits
;for degree = 2,4 do begin
        wave_fit = svdfit(fitpix/(npix/2.)-1, shift_wave, degree,/double,/legendre, yfit=pxx, measure_errors = shifterr, chisq=chisq,sigma=sigma)  
;endfor
;print,chisq

; evaluate new 1-d fit
        pixels = indgen(npix2)
        new_wave =  polyleg(pixels/(npix2/2.) - 1, wave_fit)
        new_wave2 =  polyleg(pixels/(npix2/2.) - 1, wave_fit2)

        sigma=djsig(pxx-shift)

        pixels = indgen(npix)
        new_wave_grid =  polyleg(pixels/(npix/2.) - 1, wave_fit)

    endif else begin
        message,  'WARNING: No useful bspline found. Wavelength solution not tweaked.', /info
        errflag =  1
        return, replicate(0d, n_elements(science.lambda))
    endelse

    return, new_wave2
end
