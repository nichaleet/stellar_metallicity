pro p_na, vi, na, p_rgb, p_dwf
;
; M31 RGB + MW dwarf: delta EW(Na) parameters (zone 1)
;
    na0_1d = 1.16867
    na0_1g = 0.60
    vi_12 = 1.80
    dna0_1 = 0.09
    dna_sig_1 = 1.40
    dna_sig_1_sq = dna_sig_1 * dna_sig_1
;
; M31 RGB: delta EW(Na) parameters (zones 2+3)
;
    slp_23 = -0.770649
    na0_23 = 1.97717
    ;slp_23 = -1.173913
    ;na0_23 = 3.18696
    dna0_23g = 0.10
    slp_gauss_sig = -0.8
    vi_thresh = 2.3
    dna_sig_23g = 1.00
    dna_sig_23g_sq = dna_sig_23g * dna_sig_23g
;
; MW dwarf: delta EW(Na) parameters (zones 2+3)
;
    slp_2 = 3.0605
    na0_2 = -4.37023
    vi_23 = 3.40
    na0_3 = 6.03547
    dna0_23d = 0.12
    dna_sig_23d = 1.20
    dna_sig_23d_sq = dna_sig_23d * dna_sig_23d
;
; M31 RGB: (V-I) parameters
;
    vi0_gau_g = 1.90
    vi0_exp_g = 1.23
    vi_sig_gau_g = 0.40
    vi_sig_exp_g = 0.98
    vi_frac_1g = 0.38
    vi_sig_gau_g_sq = vi_sig_gau_g * vi_sig_gau_g
;
; MW dwarf: (V-I) parameters
;
    vi0_1d = 2.56
    vi0_2d = 2.62
    vi_sig_1d = 0.37
    vi_sig_2d = 1.10
    vi_frac_1d = 0.87
    vi_sig_1d_sq = vi_sig_1d * vi_sig_1d
    vi_sig_2d_sq = vi_sig_2d * vi_sig_2d
;
; Normalization:
;
    norm_rgb = 12469.13086434567413
    norm_dwf = 14566.05061033074890
;
; Probability computation for each star in input list:
;
;
; Color distribution function: M31 RGB
;
    delvi = vi - vi0_gau_g
    delvi2 = delvi * delvi
    p_rgb_vi = vi_frac_1g * exp(-delvi2/2.0/vi_sig_gau_g_sq)
    delvi = vi - vi0_exp_g
    if delvi ge 0 then p_rgb_vi += (1.0-vi_frac_1g)*exp(-delvi/vi_sig_exp_g)
;
; Color distribution function: MW dwarf
;
    delvi = vi - vi0_1d
    delvi2 = delvi * delvi
    p_dwf_vi = vi_frac_1d * exp(-delvi2/2.0/vi_sig_1d_sq)
    delvi = vi - vi0_2d
    delvi2 = delvi * delvi
    p_dwf_vi += (1.0 - vi_frac_1d) * exp(-delvi2/2.0/vi_sig_2d_sq)

    if vi lt vi_12 then begin
;
; Na distribution function: M31 RGB (zone 1)
;
        delna = na - na0_1g - dna0_1
        delna2 = delna * delna
        f = exp(-delna2/2.0/dna_sig_1_sq) / dna_sig_1
        p_rgb = p_rgb_vi * f
;
; Na distribution function: MW dwarf (zone 1)
;
        delna = na - na0_1d - dna0_1
        delna2 = delna * delna
        f = exp(-delna2/2.0/dna_sig_1_sq) / dna_sig_1
        p_dwf = p_dwf_vi * f
    endif
;
; Na distribution function: M31 RGB (zones 2+3)
;
    if vi ge vi_12 then begin
        delta_vi = vi - vi_thresh
        if delta_vi ge 0.0 then sigma = dna_sig_23g
        if delta_vi lt 0 then sigma = dna_sig_23g + slp_gauss_sig*delta_vi
        delna = na - (vi*slp_23 + na0_23) - dna0_23g
        delna2 = delna * delna
        f = exp(-delna2/2.0/sigma/sigma) / sigma
        p_rgb = p_rgb_vi * f
;
; Na distribution function: MW dwarf (zone 2)
;
        if vi lt vi_23 then delna = na - (vi*slp_2 + na0_2) - dna0_23d
;
; Na distribution function: MW dwarf (zone 3)
;
        if vi ge vi_23 then delna = na - na0_3 - dna0_23d
        delna2 = delna * delna
        f = exp(-delna2/2.0/dna_sig_23d_sq) / dna_sig_23d
        p_dwf = p_dwf_vi * f
    endif
    p_rgb /= norm_rgb
    p_dwf /= norm_dwf

    if na gt 50.0 or na lt -90 then begin
        p_rgb = 1.0
        p_dwf = 1.0
    endif
end
