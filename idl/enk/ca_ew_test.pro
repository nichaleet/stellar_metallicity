pro gauss, x, a, f, pder
    u = (x - a[0])/a[1]
    f = a[2]*exp(-0.5*u^2) + a[3]
    pder = [[u*(f-a[3])/a[1]], [u^2.*(f-a[3])/a[1]], [(f-a[3])/a[2]], [x*0.0+1.0]]
    return
end

pro moffat, x, a, f, pder
    u = (x - a[0])/a[1]
    f = a[2]/(u^2 + 1)^a[4] + a[3]
    if n_params() ge 4 then pder = [[2.*u*a[2]*a[4]/(a[1]*(u^2 + 1)^(a[4]+1.0))], [2.*u^2.0*a[2]*a[4]/(a[1]*(u^2 + 1)^(a[4]+1.0))], [(f-a[3])/a[2]], [x*0.0+1.0], [-1.0*a[2]*alog(u^2 + 1)/(u^2 + 1)^a[4]]]
    return
end


;pro moffat2, x, a, f, pder
;    u = (x - a[0])/a[1]
;    f = a[2]/(u^2 + 1) + 1.0
;    if n_params() ge 4 then pder = [[2.*u*a[2]/(a[1]*(u^2 + 1)^(2.0))], [2.*u^2.0*a[2]/(a[1]*(u^2 + 1)^(.0))], [(f-1.0)/a[2]]]
;    return
;end

function ca_ew_test, spec1dfile, z, linepk, err=sig_ca_err_def, stop=stop
;ca_a_ew=ew_ca0, ca_b_ew=ew_ca1, ca_c_ew=ew_ca2, ca_a_ew_err=ew_ca0_err, ca_b_ew_err=ew_ca1_err, ca_c_ew_err=ew_ca2_err, cont=cont, cont_err=cont_err
    spec1d = readspec(spec1dfile)
    lam_range = [8450.0, 8700.0]
    lam_ref = [8498.43, 8542.32, 8662.44]
    frac = [0.871231, 0.397107, 0.737596]
    lambda = spec1d.lambda/(1d + z)
    lam_win = 18.0
    ;dlam = 0.33
    hwin = 0.5*lam_win
    lam_lo = lam_ref[1] + lam_win
    lam_hi = lam_ref[2] - lam_win
    lam_mid = 0.5*(lam_lo + lam_hi)
    ;N_lam_win = uint((lam_hi - lam_lo)/lam_win) + 1
    lam_cont = [lam_lo, lam_hi]

    wcont = where(lambda ge lam_cont[0] and lambda le lam_cont[1], ccont)
    ;wca = where(lambda ge lam_range[0] and lambda le lam_range[1], cca)
    wca0 = where(lambda ge (lam_ref[0] - hwin) and lambda le (lam_ref[0] + hwin), cca0)
    wca1 = where(lambda ge (lam_ref[1] - hwin) and lambda le (lam_ref[1] + hwin), cca1)
    wca2 = where(lambda ge (lam_ref[2] - hwin) and lambda le (lam_ref[2] + hwin), cca2)

;Ca A, B, C line EWs and errors
;Sig Ca and error
;Average continuum level and error


    if cca0 eq 0 or cca1 eq 0 or cca2 eq 0 then begin
        err = -999d
        return, -999d
    endif

; gaussian 1.0
; moffat 2.0
; numeric sum 3.0

    switc = 2.0


;gaussian sum
    if (switc eq 1.0) then begin
    print, 'starting gaussian sum'
    nmc = 10
    gauss_ew_ca0 = dblarr(nmc)
    gauss_ew_ca1 = dblarr(nmc)
    gauss_ew_ca2 = dblarr(nmc)
    gauss_sig_ca = dblarr(nmc)
    gcont = dblarr(nmc)
    seed = long(systime(1))
    j = 0
    for i=0,nmc-1 do begin
        broken=0.0
        j = 0
        
        try_gaussian_again:
        
        gstatus0 = 3
        gstatus1 = 3
        gstatus2 = 3
        
            spec_gcont = spec1d.spec[wcont] + (spec1d.ivar[wcont])^(-0.5)*randomn(seed, ccont)
            gcont[i] = total(spec_gcont*spec1d.ivar[wcont]) / total(spec1d.ivar[wcont])
            base = gcont[i]
        m=0
        while (gstatus0 ne 0 and broken ne 1) do begin
        
            ;why is this 1.0
            ;base = 1.0
            
            specg0 = spec1d.spec[wca0] + (spec1d.ivar[wca0])^(-0.5)*randomn(seed, cca0)
            min_specg0 = min(specg0, min_g0_ind, /nan)
            if finite(min_specg0) ne 1 then broken = 1.0
            if broken ne 1.0 then begin
                lam_r0 = lambda[wca0]
                lam_guessg0 = lam_r0[min_g0_ind]
                A_g0 = [lam_guessg0, 0.2*(max(lambda[wca0]) - min(lambda[wca0])), base - min_specg0, 0.0]
                fit_g0 = curvefit(lambda[wca0], base - specg0, spec1d.ivar[wca0], A_g0, A_g0_err, function_name='gauss', chisq=chisq_gauss0, fita=[1,1,1,0], status=gstatus0)
            endif
            
            ;if finite(min_specg0) ne 1 then begin
            ;    gstatus0 = 0
            ;    A_g0 = [lam_guessg0, 0.2*(max(lambda[wca0]) - min(lambda[wca0])), base, 0.0]
            ;    broken = 1
            ;endif
            ;print, 'gstatus0', gstatus0
            ;spec_cont = spec1d.spec[wcont]
            m = m + 1
            if m ge 100 then broken = 1.0
            
        endwhile
        
        ;print, 'line a', gstatus0, A_g0, i
        m=0
        while (gstatus1 ne 0 and broken ne 1) do begin
        
            ;base = 1.0
            specg1 = spec1d.spec[wca1] + (spec1d.ivar[wca1])^(-0.5)*randomn(seed, cca1)
            min_specg1 = min(specg1, min_g1_ind, /nan)
            if finite(min_specg1) ne 1 then broken = 1.0
            if broken ne 1 then begin
                lam_r1 = lambda[wca1]
                lam_guessg1 = lam_r1[min_g1_ind]
                A_g1 = [lam_guessg1, 0.2*(max(lambda[wca1]) - min(lambda[wca1])), base - min_specg1, 0.0]
                fit_g1 = curvefit(lambda[wca1], base - specg1, spec1d.ivar[wca1], A_g1, A_g1_err, function_name='gauss', chisq=chisq_gauss1, fita=[1,1,1,0], status=gstatus1)
            endif
            
            ;if finite(min_specg1) ne 1 then begin
            ;    gstatus1 = 0
            ;    A_g1 = [lam_guessg1, 0.2*(max(lambda[wca1]) - min(lambda[wca1])), base, 0.0]
            ;    broken = 1
            ; endif
            print, 'gstatus1', gstatus1
            m = m + 1
            if m ge 100 then broken = 1.0
            
        endwhile
        
        ;print, 'line b', gstatus1, A_g1, i
        
        m=0
        while (gstatus2 ne 0 and broken ne 1) do begin
        
            ;base = 1.0
            specg2 = spec1d.spec[wca2] + (spec1d.ivar[wca2])^(-0.5)*randomn(seed, cca2)
            min_specg2 = min(specg2, min_g2_ind, /nan)
            if finite(min_specg2) ne 1 then broken = 1.0
            if broken ne 1 then begin
                lam_r2 = lambda[wca2]
                lam_guessg2 = lam_r2[min_g2_ind]
                A_g2 = [lam_guessg2, 0.2*(max(lambda[wca2]) - min(lambda[wca2])), base - min_specg2, 0.0]
                fit_g2 = curvefit(lambda[wca2], base - specg2, spec1d.ivar[wca2], A_g2, A_g2_err, function_name='gauss', chisq=chisq_gauss2, fita=[1,1,1,0], status=gstatus2)
            endif
            
            ;if finite(min_specg2) ne 1 then begin
            ;    gstatus2 = 0
            ;    A_g2 = [lam_guessg2, 0.2*(max(lambda[wca2]) - min(lambda[wca2])), base, 0.0]
            ;    broken = 1
            ;endif
            ;print, 'gstatus2', gstatus2
            m = m + 1
            if m ge 100 then broken = 1.0
            
        endwhile
        
        ;print, 'line c', gstatus2, A_g2, i, broken
        
        if broken eq 1.0 then begin
            gauss_sig_ca[i] = !Values.F_NAN
            print, linepk, i, 'broken'
        endif
        
        if broken ne 1.0 then begin
            A_g0[1] = abs(A_g0[1])
            gauss_ew_ca0[i] = 1.0 * A_g0[2] * A_g0[1] * sqrt(2.*!dpi) / base
        
            A_g1[1] = abs(A_g1[1])
            gauss_ew_ca1[i] = 1.0 * A_g1[2] * A_g1[1] * sqrt(2.*!dpi) / base
        
            A_g2[1] = abs(A_g2[1])
            gauss_ew_ca2[i] = 1.0 * A_g2[2] * A_g2[1] * sqrt(2.*!dpi) / base
            
            print, 'not broken a b c ', gauss_ew_ca0[i], gauss_ew_ca1[i], gauss_ew_ca2[i] 
        
        
            magic_ew = 100.0
        
            if (gauss_ew_ca0[i] lt magic_ew) and (gauss_ew_ca1[i] lt magic_ew) and (gauss_ew_ca2[i] lt magic_ew) then $
                gauss_sig_ca[i] = 0.5 * gauss_ew_ca0[i] + gauss_ew_ca1[i] + 0.6 * gauss_ew_ca2[i] $
            else $ 
            if (gauss_ew_ca0[i] gt magic_ew) and (gauss_ew_ca1[i] lt magic_ew) and (gauss_ew_ca2[i] lt magic_ew) then $
                gauss_sig_ca[i] = (gauss_ew_ca1[i] + 0.6 * gauss_ew_ca2[i]) / frac[0] $
            else $
            if (gauss_ew_ca0[i] lt magic_ew) and (gauss_ew_ca1[i] gt magic_ew) and (gauss_ew_ca2[i] lt magic_ew) then $
                gauss_sig_ca[i] = (0.5 * gauss_ew_ca0[i] + 0.6 * gauss_ew_ca2[i]) / frac[1] $
            else $
            if (gauss_ew_ca0[i] lt magic_ew) and (gauss_ew_ca1[i] lt magic_ew) and (gauss_ew_ca2[i] gt magic_ew) then $
                gauss_sig_ca[i] = (0.5 * gauss_ew_ca0[i] + gauss_ew_ca1[i]) / frac[2] 
        ;else $
        ;if (finite(gauss_ew_ca0[i]) eq 0) then $
        ;    gauss_sig_ca[i] = -9999.9 $
        ;else $
        ;if (finite(gauss_ew_ca1[i]) eq 0) then $
        ;    gauss_sig_ca[i] = -9999.9 $
        ;else $
        ;if (finite(gauss_ew_ca2[i]) eq 0) then $
        ;    gauss_sig_ca[i] = -9999.9 $
        ;else $
        ;    gauss_sig_ca[i] = -9999.9
            
        ;if (finite(gauss_sig_ca[i]) eq 0) then begin
        ;    gauss_sig_ca[i] = -9999.9
            ;print, gauss_ew_ca0[i], gauss_ew_ca1[i], gauss_ew_ca2[i], ' ew'
        ;endif
        
            
        
        endif
        
        ;if gauss_sig_ca[i] eq 0 and i eq 0 and broken eq 0 and j gt 50 then stop
        
        if (gauss_sig_ca[i] le 0 or gauss_sig_ca[i] gt 1000.0) and (broken ne 1) and finite(gauss_sig_ca[i]) eq 1 then begin
            j = j + 1
            print, 'trying again because', gauss_sig_ca[i], i, broken, j
            if j le 1000 then goto, try_gaussian_again else begin 
                print, 'gave up ', gauss_sig_ca[i], ' line ', linepk, ' mc it ', i
                gauss_sig_ca[i] = !Values.F_NAN
                print, 'beginning next mc it'
                endelse
             
        endif
    
    endfor
    
    gauss_sig_ca_err = stddev(gauss_sig_ca, /nan)
    gauss_sig_ca = mean(gauss_sig_ca, /nan)
    
    gauss_ew_ca0_err = stddev(gauss_ew_ca0, /nan)
    gauss_ew_ca0 = mean(gauss_ew_ca0, /nan)
    gauss_ew_ca1_err = stddev(gauss_ew_ca1, /nan)
    gauss_ew_ca1 = mean(gauss_ew_ca1, /nan)
    gauss_ew_ca2_err = stddev(gauss_ew_ca2, /nan)
    gauss_ew_ca2 = mean(gauss_ew_ca2, /nan)
    
    gcont_err = stddev(gcont, /nan)
    gcont = mean(gcont, /nan)
    
    print, gauss_sig_ca, gauss_sig_ca_err, ' done'
    
    endif


;moffat sum
    if (switc eq 2.0) then begin
    print, 'starting moffat sum'
    nmc = 10
    moffat_ew_ca0 = dblarr(nmc)
    moffat_ew_ca1 = dblarr(nmc)
    moffat_ew_ca2 = dblarr(nmc)
    moffat_sig_ca = dblarr(nmc)
    mcont = dblarr(nmc)
    seed = long(systime(1))
    j = 0
    for i=0,nmc-1 do begin
        broken=0.0
        j = 0
        
        try_moffat_again:
        
        mstatus0 = 3
        mstatus1 = 3
        mstatus2 = 3
        
            spec_mcont = spec1d.spec[wcont] + (spec1d.ivar[wcont])^(-0.5)*randomn(seed, ccont)
            mcont[i] = total(spec_mcont*spec1d.ivar[wcont]) / total(spec1d.ivar[wcont])
            base = mcont[i]
        m=0
        while (mstatus0 ne 0 and broken ne 1) do begin
        
            ;why is this 1.0
            ;base = 1.0
            
            specm0 = spec1d.spec[wca0] + (spec1d.ivar[wca0])^(-0.5)*randomn(seed, cca0)
            min_specm0 = min(specm0, min_m0_ind, /nan)
            if finite(min_specm0) ne 1 then begin 
                broken = 1.0
                ;stop
            endif
            if broken ne 1.0 then begin
                lam_r0 = lambda[wca0]
                lam_guessm0 = lam_r0[min_m0_ind]
                lambda_med_m0 = median(lambda[wca0])
                ;print, lambda_med_m0, ' lambdamed'
                A_m0 = [lam_guessm0-lambda_med_m0, 0.1*(max(lambda[wca0]) - min(lambda[wca0])), (base-min_specm0)*0.8, 0.0, 1.0]
                ;print, A_m0
                fit_m0 = curvefit(lambda[wca0], base - specm0, spec1d.ivar[wca0], A_m0, A_m0_err, function_name='moffat', chisq=chisq_moffat0, fita=[0,1,1,1,0], status=mstatus0)
            endif
            
            ;if finite(min_specm0) ne 1 then begin
            ;    mstatus0 = 0
            ;    A_m0 = [lam_guessm0, 0.2*(max(lambda[wca0]) - min(lambda[wca0])), base, 0.0]
            ;    broken = 1
            ;endif
            ;print, 'mstatus0', mstatus0
            ;spec_cont = spec1d.spec[wcont]
            m = m + 1
            if m ge 100 then broken = 1.0
            
        endwhile
        
        ;print, 'line a', mstatus0, A_m0, i
        m=0
        while (mstatus1 ne 0 and broken ne 1) do begin
        
            ;base = 1.0
            specm1 = spec1d.spec[wca1] + (spec1d.ivar[wca1])^(-0.5)*randomn(seed, cca1)
            min_specm1 = min(specm1, min_m1_ind, /nan)
            if finite(min_specm1) ne 1 then broken = 1.0
            if broken ne 1 then begin
                lam_r1 = lambda[wca1]
                lam_guessm1 = lam_r1[min_m1_ind]
                lambda_med_m1 = median(lambda[wca1])
                A_m1 = [lam_guessm1-lambda_med_m1, 0.1*(max(lambda[wca1]) - min(lambda[wca1])), (base-min_specm1)*0.8, 0.0, 1.0]
                fit_m1 = curvefit(lambda[wca1], base - specm1, spec1d.ivar[wca1], A_m1, A_m1_err, function_name='moffat', chisq=chisq_moffat1, fita=[0,1,1,1,0], status=mstatus1)
            endif
            
            ;if finite(min_specm1) ne 1 then begin
            ;    mstatus1 = 0
            ;    A_m1 = [lam_guessm1, 0.2*(max(lambda[wca1]) - min(lambda[wca1])), base, 0.0]
            ;    broken = 1
            ; endif
            ;print, 'mstatus1', mstatus1
            m = m + 1
            if m ge 100 then broken = 1.0
            
        endwhile
        
        ;print, 'line b', mstatus1, A_m1, i
        
        m=0
        while (mstatus2 ne 0 and broken ne 1) do begin
        
            ;base = 1.0
            specm2 = spec1d.spec[wca2] + (spec1d.ivar[wca2])^(-0.5)*randomn(seed, cca2)
            min_specm2 = min(specm2, min_m2_ind, /nan)
            if finite(min_specm2) ne 1 then broken = 1.0
            if broken ne 1 then begin
                lam_r2 = lambda[wca2]
                lam_guessm2 = lam_r2[min_m2_ind]
                lambda_med_m2 = median(lambda[wca2])
                A_m2 = [lam_guessm2-lambda_med_m2, 0.1*(max(lambda[wca2]) - min(lambda[wca2])), (base-min_specm2)*0.8, 0.0, 1.0]
                fit_m2 = curvefit(lambda[wca2], base - specm2, spec1d.ivar[wca2], A_m2, A_m2_err, function_name='moffat', chisq=chisq_moffat2, fita=[0,1,1,1,0], status=mstatus2)
            endif
            
            ;if finite(min_specm2) ne 1 then begin
            ;    mstatus2 = 0
            ;    A_m2 = [lam_guessm2, 0.2*(max(lambda[wca2]) - min(lambda[wca2])), base, 0.0]
            ;    broken = 1
            ;endif
            ;print, 'mstatus2', mstatus2
            m = m + 1
            if m ge 100 then broken = 1.0
            
        endwhile
        
        ;print, 'line c', mstatus2, A_m2, i, broken
        
        if broken eq 1.0 then begin
            moffat_sig_ca[i] = !Values.F_NAN
            print, linepk, i, 'broken'
        endif
        
        if broken ne 1.0 then begin
        A_m0[1] = abs(A_m0[1])
        moffat_ew_ca0[i] = 1.0 * A_m0[2] * A_m0[1] * !dpi / base
        
        A_m1[1] = abs(A_m1[1])
        moffat_ew_ca1[i] = 1.0 * A_m1[2] * A_m1[1] * !dpi / base
        
        A_m2[1] = abs(A_m2[1])
        moffat_ew_ca2[i] = 1.0 * A_m2[2] * A_m2[1] * !dpi / base
        
        blah = A_m0[2] * A_m0[1] * 3.1415927 / base
            
            print, A_m0
            print, A_m0[2], A_m0[1], base
            print, blah
            ;print, 'not broken a b c ', moffat_ew_ca0[i], moffat_ew_ca1[i], moffat_ew_ca2[i] 
        
        
            magic_ew = 100.0
        
            if (moffat_ew_ca0[i] lt magic_ew) and (moffat_ew_ca1[i] lt magic_ew) and (moffat_ew_ca2[i] lt magic_ew) then $
                moffat_sig_ca[i] = 0.5 * moffat_ew_ca0[i] + moffat_ew_ca1[i] + 0.6 * moffat_ew_ca2[i] $
            else $ 
            if (moffat_ew_ca0[i] gt magic_ew) and (moffat_ew_ca1[i] lt magic_ew) and (moffat_ew_ca2[i] lt magic_ew) then $
                moffat_sig_ca[i] = (moffat_ew_ca1[i] + 0.6 * moffat_ew_ca2[i]) / frac[0] $
            else $
            if (moffat_ew_ca0[i] lt magic_ew) and (moffat_ew_ca1[i] gt magic_ew) and (moffat_ew_ca2[i] lt magic_ew) then $
                moffat_sig_ca[i] = (0.5 * moffat_ew_ca0[i] + 0.6 * moffat_ew_ca2[i]) / frac[1] $
            else $
            if (moffat_ew_ca0[i] lt magic_ew) and (moffat_ew_ca1[i] lt magic_ew) and (moffat_ew_ca2[i] gt magic_ew) then $
                moffat_sig_ca[i] = (0.5 * moffat_ew_ca0[i] + moffat_ew_ca1[i]) / frac[2] 
        ;else $
        ;if (finite(moffat_ew_ca0[i]) eq 0) then $
        ;    moffat_sig_ca[i] = -9999.9 $
        ;else $
        ;if (finite(moffat_ew_ca1[i]) eq 0) then $
        ;    moffat_sig_ca[i] = -9999.9 $
        ;else $
        ;if (finite(moffat_ew_ca2[i]) eq 0) then $
        ;    moffat_sig_ca[i] = -9999.9 $
        ;else $
        ;    moffat_sig_ca[i] = -9999.9
            
        ;if (finite(moffat_sig_ca[i]) eq 0) then begin
        ;    moffat_sig_ca[i] = -9999.9
            ;print, moffat_ew_ca0[i], moffat_ew_ca1[i], moffat_ew_ca2[i], ' ew'
        ;endif
        
            
        
        endif
        
        ;if moffat_sig_ca[i] eq 0 and i eq 0 and broken eq 0 and j gt 50 then stop
        
        if (moffat_sig_ca[i] le 0 or moffat_sig_ca[i] gt 1000.0) and (broken ne 1) and finite(moffat_sig_ca[i]) eq 1 then begin
            j = j + 1
            print, 'trying again because', moffat_sig_ca[i], i, broken, j
            if j le 1000 then goto, try_moffat_again else begin 
                print, 'gave up ', moffat_sig_ca[i], ' line ', linepk, ' mc it ', i
                moffat_sig_ca[i] = !Values.F_NAN
                print, 'beginning next mc it'
                endelse
             
        endif
    
    endfor
    
    moffat_sig_ca_err = stddev(moffat_sig_ca, /nan)
    moffat_sig_ca = mean(moffat_sig_ca, /nan)
    
    moffat_ew_ca0_err = stddev(moffat_ew_ca0, /nan)
    moffat_ew_ca0 = mean(moffat_ew_ca0, /nan)
    moffat_ew_ca1_err = stddev(moffat_ew_ca1, /nan)
    moffat_ew_ca1 = mean(moffat_ew_ca1, /nan)
    moffat_ew_ca2_err = stddev(moffat_ew_ca2, /nan)
    moffat_ew_ca2 = mean(moffat_ew_ca2, /nan)
    
    mcont_err = stddev(mcont, /nan)
    mcont = mean(mcont, /nan)
    
    print, moffat_sig_ca, moffat_sig_ca_err, 'done'
    
    endif

    if (switc eq 3.0) then begin
;ivar-weighted numeric sum
    nmc = 1000
    ew_ca0 = dblarr(nmc)
    ew_ca1 = dblarr(nmc)
    ew_ca2 = dblarr(nmc)
    sig_ca = dblarr(nmc)
    cont = dblarr(nmc)
    seed = long(systime(1))
    for i=0,nmc-1 do begin
    
        spec0 = spec1d.spec[wca0] + (spec1d.ivar[wca0])^(-0.5)*randomn(seed, cca0)
        spec1 = spec1d.spec[wca1] + (spec1d.ivar[wca1])^(-0.5)*randomn(seed, cca1)
        spec2 = spec1d.spec[wca2] + (spec1d.ivar[wca2])^(-0.5)*randomn(seed, cca2)
        
        ;spec0 = spec1d.spec[wca0]
        ;spec1 = spec1d.spec[wca1]
        ;spec2 = spec1d.spec[wca2]
        
        spec_cont = spec1d.spec[wcont] + (spec1d.ivar[wcont])^(-0.5)*randomn(seed, ccont)
        cont[i] = total(spec_cont*spec1d.ivar[wcont]) / total(spec1d.ivar[wcont])
        
        ;spec_cont = spec1d.spec[wcont]
        
        wcas0 = [wca0, wca0[cca0-1]+1]
        dlam0 = shift(lambda[wcas0],-1) - lambda[wcas0]
        dlam0 = dlam0[0:cca0-1]
        
        wcas1 = [wca1, wca1[cca1-1]+1]
        dlam1 = shift(lambda[wcas1],-1) - lambda[wcas1]
        dlam1 = dlam1[0:cca1-1]

        wcas2 = [wca2, wca2[cca2-1]+1]
        dlam2 = shift(lambda[wcas2],-1) - lambda[wcas2]
        dlam2 = dlam2[0:cca2-1]
        
        ew_ca0[i] = cca0 * total( ( (cont - spec0)*spec1d.ivar[wca0]*dlam0 ) / cont ) / total(spec1d.ivar[wca0])
        ew_ca1[i] = cca1 * total( ( (cont - spec1)*spec1d.ivar[wca1]*dlam1 ) / cont ) / total(spec1d.ivar[wca1])
        ew_ca2[i] = cca2 * total( ( (cont - spec2)*spec1d.ivar[wca2]*dlam2 ) / cont ) / total(spec1d.ivar[wca2])
        
        ;ew_ca0[i] = total( ( (cont - spec0)*dlam0 ) / cont )
        ;ew_ca1[i] = total( ( (cont - spec1)*dlam1 ) / cont )
        ;ew_ca2[i] = total( ( (cont - spec2)*dlam2 ) / cont )
        
        magic_ew = 100.0
        
        if (ew_ca0[i] lt magic_ew) and (ew_ca1[i] lt magic_ew) and (ew_ca2[i] lt magic_ew) then $
            sig_ca[i] = 0.5 * ew_ca0[i] + ew_ca1[i] + 0.6 * ew_ca2[i] $
        else $ 
        if (ew_ca0[i] gt magic_ew) and (ew_ca1[i] lt magic_ew) and (ew_ca2[i] lt magic_ew) then $
            sig_ca[i] = (ew_ca1[i] + 0.6 * ew_ca2[i]) / frac[0] $
        else $
        if (ew_ca0[i] lt magic_ew) and (ew_ca1[i] gt magic_ew) and (ew_ca2[i] lt magic_ew) then $
            sig_ca[i] = (0.5 * ew_ca0[i] + 0.6 * ew_ca2[i]) / frac[1] $
        else $
        if (ew_ca0[i] lt magic_ew) and (ew_ca1[i] lt magic_ew) and (ew_ca2[i] gt magic_ew) then $
            sig_ca[i] = (0.5 * ew_ca0[i] + ew_ca1[i]) / frac[2] 
        ;else $
        ;if (finite(ew_ca0[i]) eq 0) then $
        ;    sig_ca[i] = -9999.9 $
        ;else $
        ;if (finite(ew_ca1[i]) eq 0) then $
        ;    sig_ca[i] = -9999.9 $
        ;else $
        ;if (finite(ew_ca2[i]) eq 0) then $
        ;    sig_ca[i] = -9999.9 $
        ;else $
        ;    sig_ca[i] = -9999.9
            
        ;if (finite(sig_ca[i]) eq 0) then begin
        ;    sig_ca[i] = -9999.9
            ;print, ew_ca0[i], ew_ca1[i], ew_ca2[i], ' ew'
        ;endif
        
        ;print, sig_ca[i]
                    
    endfor
    
    sig_ca_err = stddev(sig_ca, /nan)
    sig_ca = mean(sig_ca, /nan)
    
    ew_ca0_err = stddev(ew_ca0, /nan)
    ew_ca0 = mean(ew_ca0, /nan)
    ew_ca1_err = stddev(ew_ca1, /nan)
    ew_ca1 = mean(ew_ca1, /nan)
    ew_ca2_err = stddev(ew_ca2, /nan)
    ew_ca2 = mean(ew_ca2, /nan)
    
    cont_err = stddev(cont, /nan)
    cont = mean(cont, /nan)
    
    endif

    sig_ca_err_def = moffat_sig_ca_err
    print, 'REALLY done, back to consolidate'
    return, moffat_sig_ca
end