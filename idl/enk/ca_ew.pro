function ca_ew, spec1dfile, z, err=sig_ca_err, stop=stop
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
        sig_ca_err = -999d
        return, -999d
    endif

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

    return, sig_ca
end
