function Pmcmc, vobs, sobs, v1, s1
    N = n_elements(vobs)
    P = -0.5*total(alog(sobs^2 + s1^2)) - 0.5*total((vobs-v1)^2 / (sobs^2 + s1^2)) - N*alog(2.*!PI)/2.
    return, P
end


function ek_mcmc_sigma, vobs, sobs,v0, s0, m, m_low, m_high, ntrial=ntrial, var=var
    ; SET UP INITIAL RANGE OVER WHICH TO SEARCH
    if ~keyword_set(var) then var=2.5

    v1 = v0
    s1 = s0

    vavg = v1
    j=0L

    if ~keyword_set(ntrial) then ntrial = long(1e6)
    s_good = fltarr(ntrial)
    v_good = fltarr(ntrial)

    r1=randomn(seed,ntrial)
    r2=randomn(seed+1,ntrial)
    r3=randomu(seed+2,ntrial)

    P1 = Pmcmc(vobs, sobs, v1, s1)

    for ii=0,ntrial-2 do begin
        ; Sample next step
        s2 = abs(s1 + r1[ii]*var)
        v2 = v1 + r2[ii]*var

        ; Calculate second prob
        P2 = Pmcmc(vobs, sobs, v2, s2)

        ; Accept or reject parameters
        P = exp(P2-P1)
        if (P ge 1) or (r3[ii] lt P) then begin
            v1 = v2
            s1 = s2
            P1 = P2
            v_good[j] = v1
            s_good[j] = s1
            j++
        endif
    endfor   
    if (j le 100) then begin
        message,'not enough accepts, increase loop'
    endif
    v_good = v_good[0:j-1]
    s_good = s_good[0:j-1]

; ERRORS on SIGMA
    cgHistoplot,s_good,/oprob,probability=p,locations=sval,bin=0.01
    n=n_elements(p)
    bin = Value_Locate(p, 0.8415)
    c1s=bin+1
    bin = Value_Locate(p, 0.1585)
    c2s=bin

    bin = Value_Locate(p, 0.5)
    if Abs(p[bin] - 0.5) GT Abs(p[bin+1] - 0.5) THEN $
        c3s = bin+1 ELSE c3s = bin

    err = (sval[c1s] -sval[c2s])/2.
;print,'sigma =     ',sval[c3s], ' +/- ',err

    m=sval[c3s]
    m_low = sval[c3s]-sval[c2s]
    m_high = sval[c1s]-sval[c3s]
;print,'Upper/Lower Limit = ',sval[c3s],' + ',sval[c1s]-sval[c3s],' - ',sval[c3s]-sval[c2s]
;print


; ERRORS on VELOCITY
    cgHistoplot,v_good,/oprob,probability=p,locations=vval,bin=0.05
    n=n_elements(p)
    bin = Value_Locate(p, 0.8415)
    c1v=bin+1
    bin = Value_Locate(p, 0.1585)
    c2v=bin

    bin = Value_Locate(p, 0.5)
    if Abs(p[bin] - 0.5) GT Abs(p[bin+1] - 0.5) THEN $
        c3v = bin+1 ELSE c3v = bin

    err = (vval[c1v] -vval[c2v])/2.
;print,'vel =      ',vval[c3v], ' +/- ',err



    return, {v0:vval[c3v], v0err:(vval[c1v]-vval[c2v])/2., v0errlower:vval[c3v]-vval[c2v], v0errupper:vval[c1v]-vval[c3v], sigmav:sval[c3s], sigmaverr:(sval[c1s]-sval[c2s])/2., sigmaverrlower:sval[c3s]-sval[c2s], sigmaverrupper:sval[c1s]-sval[c3s], v0arr:v_good, sigmavarr:s_good, ntrials:j}

end
