pro p_cmd, x, y, p_rgb, p_dwf
;
; M31 RGB y parameters:
;
    y0_g = 0.510
    y_sig_g = 0.310
    y_sig_g_sq = y_sig_g * y_sig_g
;
; M31 RGB: x paramters: 
;
    gauss_mn_a_g = 0.189451
    gauss_mn_b_g = -0.1572680
    gauss_mn_c_g = 0.308648
    gauss_sig_a_g= 0.220984
    gauss_sig_b_g= -0.2546160
    gauss_sig_c_g= 0.222628
;
; MW dwarf: y paramters 
;
    y0_d = 0.810
    y_sig_d = 0.330
    y_sig_d_sq = y_sig_d * y_sig_d
;
; MW dwarf: x parameters:
;
    gauss_mn_a_d =-0.158741
    gauss_mn_b_d =0.0305537
    gauss_mn_c_d =0.478445
    gauss_sig_a_d=-0.101776
    gauss_sig_b_d=0.0295984
    gauss_sig_c_d=0.224484
    y3=1.06804
    x0_d_y3 = 0.330
    x_sig_d_y3 = 0.140
    slp_mn_d = 2.0*gauss_mn_a_d*y3 + gauss_mn_b_d
    slp_sig_d = 2.0*gauss_sig_a_d*y3 + gauss_sig_b_d
;
; Normalization:
;
    norm_rgb = 19643.9947893616
    norm_dwf = 20734.5147046192
;
; Probability computation for each star:
;
;
; y distribution function: M31 RGB
;
    y2=y*y
    dely = y - y0_g
    dely2 = dely * dely
    p_rgb_y = exp(-dely2/2.0/y_sig_g_sq)
;
; y distribution function: MW dwarf
;
    dely = y - y0_d
    dely2 = dely * dely
    p_dwf_y = exp(-dely2/2.0/y_sig_d_sq)
;
; Compute x Gaussian parameters for giants
;
    x0_g=gauss_mn_a_g*y2+gauss_mn_b_g*y+gauss_mn_c_g
    x_sig_g=gauss_sig_a_g*y2+gauss_sig_b_g*y+gauss_sig_c_g
    x_sig_g_sq=x_sig_g*x_sig_g
;
; Compute x Gaussian parameters for dwarfs
;
    x0_d=gauss_mn_a_d*y2+gauss_mn_b_d*y+gauss_mn_c_d
    x_sig_d=gauss_sig_a_d*y2+gauss_sig_b_d*y+gauss_sig_c_d
;
; Avoid quadratic extrapolation: linear extrapolation for gaussian mean
;          constant value extrapolation for gaussian sigma
;
    if y gt y3 then begin
        delta_y = y - y3
        x0_d = x0_d_y3 + delta_y * slp_mn_d
        ;x_sig_d = x_sig_d_y3 + delta_y * slp_sig_d
        ;x0_d = x0_d_y3
        x_sig_d = x_sig_d_y3
    endif
    x_sig_d_sq=x_sig_d*x_sig_d
;
; x distribution function: M31 RGB 
;
    delx = x - x0_g
    delx2 = delx * delx
    f = exp(-delx2/2.0/x_sig_g_sq) / x_sig_g
    p_rgb = p_rgb_y * f
;
; x distribution function: MW dwarf 
;
    delx = x - x0_d
    delx2 = delx * delx
    f = exp(-delx2/2.0/x_sig_d_sq) / x_sig_d
    p_dwf = p_dwf_y * f
    
    p_rgb /= norm_rgb
    p_dwf /= norm_dwf

    ;if x ge 10 or x le -10 or y ge 10 or y le -10
    if x ge 9.99 or x le -10 or y ge 9.99 or y le -10 then begin
        p_rgb = 1.00
        p_dwf = 1.00
    endif
end
