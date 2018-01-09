pro p_ddo, ddo, p_rgb, p_dwf
;
; RGB: Single exponential parameters / Normalization
;
    ddo0 = 0.25
    ddo_hi = 1.0
    norm_rgb = 0.24542
;
; Dwarf: Double exponential parameters / Normalization
;
    ddo0_1 = 0.04
    ddo0_2 = 0.60
    frac_1 = 0.88
    norm_dwf = 0.0935973
    
;
; Probability computation for each star:
;
    if ddo gt 1.0 and ddo le 9 then ddo = 1.0
    if ddo le 0.0 then ddo = 0.0

    p_rgb = exp((ddo-ddo_hi)/ddo0)/norm_rgb
    f_1 = frac_1*exp(-ddo/ddo0_1)
    f_2 = (1.0 - frac_1)*exp(-ddo/ddo0_2)
    p_dwf = (f_1 + f_2)/norm_dwf

    if ddo gt 9 then begin           
        p_rgb = 1.0
        p_dwf = 1.0
    endif
end
