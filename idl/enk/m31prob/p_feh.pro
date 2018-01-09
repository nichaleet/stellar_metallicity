pro p_feh, feh_phot, feh_spec, p_rgb, p_dwf
    max_exp_arg = 725.0
;
; M31 RGB: Fe/H_phot (y) parameters:
;
    y0_g_1 = -0.30
    y_sig_g_1 = 0.10
    y_sig_g_1_sq = y_sig_g_1 * y_sig_g_1
    y0_g_2 = -0.60
    y_sig_g_2 = 0.40
    y_sig_g_2_sq = y_sig_g_2 * y_sig_g_2
    frac_g= 0.7
;
; M31 RGB: Fe/H_spec (x) parameters: 
;
    gauss_mn_a_g = 0.741283
    gauss_mn_b_g = 2.184190
    gauss_mn_c_g = 0.324958
    gauss_sig_a_g= 0.120906
    gauss_sig_b_g= 0.134015
    gauss_sig_c_g= 0.748239
;
; MW dwarf: Fe/H_phot (y) parameters:
;
    y0_d_1 = -0.28
    y_sig_d_1 = 0.32
    y0_d_2 = -0.70
    y_sig_d_2 = 0.40
    y0_d_3 = -2.10
    y_sig_d_3 = 0.28
    frac_d_1 = 0.90 
    frac_d_2 = 0.02
    y_exp = 0.2

    y_sig_d_1_sq = y_sig_d_1 * y_sig_d_1
    y_sig_d_2_sq = y_sig_d_2 * y_sig_d_2
    y_sig_d_3_sq = y_sig_d_3 * y_sig_d_3
;
; MW dwarf: Fe/H_spec (x) parameters:
;
    gauss_mn_a_d  =  0.241325
    gauss_mn_b_d  =  0.634522
    gauss_mn_c_d  = -1.559290
    gauss_sig_a_d =  0.290000
    gauss_sig_b_d =  0.599769
    gauss_sig_c_d =  0.702870
;
; Normalization:
;
    norm_rgb = 11932.8
    norm_dwf = 19218.1
    
    
;  Fe/H_spec:
    x = feh_spec
;  Fe/H_phot:
    y = feh_phot

    y2=y*y
;
; Fe/H_phot distribution function: M31 RGB
;
    dely_1 = y - y0_g_1
    dely_2 = y - y0_g_2
    dely2_1 = dely_1 * dely_1
    dely2_2 = dely_2 * dely_2
    t = dely2_1/2.0/y_sig_g_1_sq
    if t ge max_exp_arg then t = max_exp_arg
    p_rgb_y = frac_g*exp(-1.0*t)
    t = dely2_2/2.0/y_sig_g_2_sq
    if t ge max_exp_arg then t = max_exp_arg
    p_rgb_y += (1.0 - frac_g)*exp(-1.0*t)
;
; Fe/H_phot distribution function: MW dwarf
;
    dely_1 = y - y0_d_1
    dely_2 = y - y0_d_2
    dely_3 = y - y0_d_3
    dely2_1 = dely_1 * dely_1
    dely2_2 = dely_2 * dely_2
    dely2_3 = dely_3 * dely_3
    t = dely2_1/2.0/y_sig_d_1_sq
    if t ge max_exp_arg then t = max_exp_arg
    p_dwf_y = frac_d_1*exp(-1.0*t)
    t = dely2_2/2.0/y_sig_d_2_sq
    if t ge max_exp_arg then t = max_exp_arg
    p_dwf_y += (1.0 - frac_d_1)*exp(-1.0*t)
    t = dely2_3/2.0/y_sig_d_3_sq
    if t ge max_exp_arg then t = max_exp_arg
    p_dwf_y += frac_d_2*exp(-1.0*t)
    t = y/y_exp
    if t ge max_exp_arg then t = max_exp_arg
    if y ge 0 then p_dwf_y *= exp(-1.0*t)
;
;
; Compute Fe/H_spec Gaussian parameters: M31 RGB
;
    x0_g=gauss_mn_a_g*y2+gauss_mn_b_g*y+gauss_mn_c_g
    x_sig_g=gauss_sig_a_g*y2+gauss_sig_b_g*y+gauss_sig_c_g
    x_sig_g_sq=x_sig_g*x_sig_g
;
; Compute Fe/H_spec Gaussian parameters: MW dwarf
;
    x0_d=gauss_mn_a_d*y2+gauss_mn_b_d*y+gauss_mn_c_d
    x_sig_d=gauss_sig_a_d*y2+gauss_sig_b_d*y+gauss_sig_c_d
    x_sig_d_sq=x_sig_d*x_sig_d

;
; Fe/H_spec distribution function: M31 RGB 
;
    delx = x - x0_g
    delx2 = delx * delx
    t = delx2/2.0/x_sig_g_sq
    if t ge max_exp_arg then t = max_exp_arg
    f = exp(-1.0*t) / x_sig_g
    p_rgb = p_rgb_y * f
;
; Fe/H_spec distribution function: MW dwarf 
;
    delx = x - x0_d
    delx2 = delx * delx
    t = delx2/2.0/x_sig_d_sq
    if t ge max_exp_arg then t = max_exp_arg
    f = exp(-1.0*t) / x_sig_d
    p_dwf = p_dwf_y * f
    
    p_rgb /= norm_rgb
    p_dwf /= norm_dwf

    if x le -10 || x ge 10 || y le -10 || y ge 10 then begin
    ;if (x le -5 || x ge 3 || y le -5 || y ge 3)
        p_rgb = 1.0
        p_dwf = 1.0
    endif
end
