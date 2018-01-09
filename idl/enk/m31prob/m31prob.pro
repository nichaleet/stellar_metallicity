function m31prob, vel, gprob, vi, na, x, y, feh_phot, feh_spec, class=class, savflag=savflag
    if ~keyword_set(savflag) then savflag = 'm31'

    ndiag = 5
    n = n_elements(vel)
    
    log10 = alog(10.0)
    Lcap = 5.0
;
; SELECT diagnostics
;
    L = dblarr(n)
    class = intarr(n)
    for j=0,n-1 do begin
        flg = intarr(ndiag)
        flg[0:4] = 1
        if ndiag gt 5 then flg[5:ndiag-1] = 0
;
; SET ERROR LEVELS in the probs
;
        epd = dblarr(ndiag)
        ;epd[0] = 0.00001
        epd[0] = 0.0003
        epd[1] = 0.2
        ;epd[2] = 0.0000003
        epd[2] = 0.000003
        ;epd[3] = 0.000003
        epd[3] = 0.00003
        ;epd[4] = 0.0000005
        epd[4] = 0.000005
        if ndiag gt 5 then epd[5:ndiag-1] = 1.0

        epg = dblarr(ndiag)
        ;epg[0] = 0.00001
        epg[0] = 0.0003
        epg[1] = 0.05
        ;epg[2] = 0.0000003
        epg[2] = 0.000003
        ;epg[3] = 0.000002
        epg[3] = 0.00002
        ;epg[4] = 0.0000005
        epg[4] = 0.000005
        if ndiag gt 5 then epg[5:ndiag-1] = 1.0
;
; INITIALIZE probs
;
        pd = dblarr(ndiag)+1.0
        pg = dblarr(ndiag)+1.0

        ;if keyword_set(novel) then flg[0] = 0 else flg[0] = 1
        ;if flg[0] eq 1 then begin
            p_vel, vel[j], pg0, pd0, savflag=savflag
            pg[0] = pg0  &  pd[0] = pd0
        ;endif

        if abs(gprob[j]) ge 9.99 or ~finite(gprob[j]) then flg[1] = 0 else flg[1] = 1
        if flg[1] eq 1 then begin
            p_ddo, gprob[j], pg1, pd1
            pg[1] = pg1  &  pd[1] = pd1
        endif

        if na[j] ge 99 or na[j] le -9 or ~finite(na[j]) then flg[2] = 0 else flg[2] = 1
        if flg[2] eq 1 then begin
            p_na, vi[j], na[j], pg2, pd2
            pg[2] = pg2  &  pd[2] = pd2
        endif

        if ~finite(x[j]) or ~finite(y[j]) then flg[3] = 0 else flg[3] = 1
        if flg[3] eq 1 then begin
            p_cmd, x[j], y[j], pg3, pd3
            pg[3] = pg3  &  pd[3] = pd3
        endif

        if abs(feh_spec[j]) gt 9.99 or ~finite(feh_phot[j]) or ~finite(feh_spec[j]) then flg[4] = 0 else flg[4] = 1
        if flg[4] eq 1 then begin
            p_feh, feh_phot[j], feh_spec[j], pg4, pd4
            pg[4] = pg4  &  pd[4] = pd4
        endif

        sn=0.0
        sd=0.0
        Lij = dblarr(ndiag)
        for i=0,ndiag-1 do begin
            w = double(flg[i])
            rg = epg[i] * epg[i] / pg[i] / pg[i]
            rd = epd[i] * epd[i] / pd[i] / pd[i]
            if i ne 0 and i ne 1 and rg ge 0.1 and rd ge 0.1 then w /= 5.0*(rg+rd)
            Lij = alog(pg[i]/pd[i]) / log10
            if Lij ge Lcap then Lij = Lcap
            if Lij le -1.0*Lcap then Lij = -1.0*Lcap
            sn += w*Lij
            sd += w
        endfor
        L[j] = sn/sd

        case 1 of
            L[j] gt 0.5 and x[j] ge 0: class[j] = 3
            L[j] gt 0.5 and x[j] ge -0.05 and x[j] lt 0: class[j] = 2
            L[j] le 0.5 and L[j] ge 0 and x[j] ge -0.05: class[j] = 1
            L[j] ge -0.5 and L[j] lt 0 and x[j] ge -0.05: class[j] = -1
            L[j] ge -0.5 and x[j] lt -0.05: class[j] = -2
            L[j] lt -0.5: class[j] = -3
            else: class[j] = 0
        endcase
    endfor

    return, L
end
