pro p_vel, vel, p_rgb, p_dwf, savflag=savflag
    if ~keyword_set(veldist) then veldist = 0
;
; Dwarf: Double Gaussian parameters:
;
    vel_mean1 = -38
    vel_mean2 =-95
    vel_sig1 = 30
    vel_sig2 = 48
    frac1 = 0.56

    vel_sig1_sq = vel_sig1*vel_sig1
    vel_sig2_sq = vel_sig2*vel_sig2
;
; RGB: Single Gaussian parameters:
;
    case savflag of
        'and1': begin
            vel_mean = -366.4
            vel_sig = 45.7
        end
        'and2': begin
            vel_mean = -199.4
            vel_sig = 16.0
        end
        'and3': begin
            vel_mean = -314.4
            vel_sig = 70.4
        end
        else: begin
            vel_mean = -300.
            vel_sig = 85.
        end

    endcase

    vel_sig2 = vel_sig*vel_sig
;
; Normalization:
;
    norm_rgb = 213.046
    norm_dwf = 95.05

    delvel1 = vel - vel_mean1
    delvel2 = vel - vel_mean2
    p_dwf = frac1*exp(-1.0*delvel1*delvel1/2.0/vel_sig1_sq)
    p_dwf += (1.0 - frac1) * exp(-1.0*delvel2*delvel2/2.0/vel_sig2_sq)
    p_dwf /= norm_dwf

    delvel = vel - vel_mean
    p_rgb = exp(-1.0*delvel*delvel/2.0/vel_sig2)
    p_rgb /= norm_rgb
end
