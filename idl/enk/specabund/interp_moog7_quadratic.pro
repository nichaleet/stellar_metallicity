function npar, array
    return, round((double(array[1])-double(array[0]))/double(array[2])) + 1
end


function makepargrid, array, npar
    return, dindgen(npar)*array[2] + array[0]
end


function interp_moog7, teff, logg, feh, alpha, mooglambda=mooglambda, fullres=fullres
    teff_grid1 = [3500., 5501., 100.]
    teff_grid2 = [5600., 8001., 200.]
    logg_grid = [0.0, 5.01, 0.5]
    feh_grid = [-5.0, 0.01, 0.1]
    alpha_grid = [-0.8, 1.21, 0.1]

    nteff1 = npar(teff_grid1)
    nteff2 = npar(teff_grid2)
    nteff = nteff1 + nteff2
    nlogg = npar(logg_grid)
    nfeh = npar(feh_grid)
    nalpha = npar(alpha_grid)

    ateff = [makepargrid(teff_grid1, nteff1), makepargrid(teff_grid2, nteff2)]
    alogg = makepargrid(logg_grid, nlogg)
    afeh = makepargrid(feh_grid, nfeh)
    aalpha = makepargrid(alpha_grid, nalpha)

    eps = 1d-7

    dteff = ateff - teff
    w0 = where(abs(dteff) lt eps, c0)
    if c0 eq 1 then begin
        dteff = 1d
        tteff = ateff[w0]
        iteff = 1
    endif else begin
        wpos = where(dteff gt 0)
        junk = min(dteff[wpos], wteff)
        if (abs(dteff[wpos[wteff]-1]) lt abs(dteff[wpos[wteff]]) and wpos[wteff]-1 gt 0) or (wpos[wteff] ge nteff-1) then begin
            dteff = reverse(dteff[wpos[wteff]-2:wpos[wteff]])
            tteff = ateff[wpos[wteff]-2:wpos[wteff]]
        endif else begin
            dteff = reverse(dteff[wpos[wteff]-1:wpos[wteff]+1])
            tteff = ateff[wpos[wteff]-1:wpos[wteff]+1]
        endelse
        iteff = 3
    endelse

    dlogg = alogg - logg
    w0 = where(abs(dlogg) lt eps, c0)
    if c0 eq 1 then begin
        dlogg = 1d
        tlogg = alogg[w0]
        ilogg = 1
    endif else begin
        wpos = where(dlogg gt 0)
        junk = min(dlogg[wpos], wlogg)
        if (abs(dlogg[wpos[wlogg]-1]) lt abs(dlogg[wpos[wlogg]]) and wpos[wlogg]-1 gt 0) or (wpos[wlogg] ge nlogg) then begin
            dlogg = reverse(dlogg[wpos[wlogg]-2:wpos[wlogg]])
            tlogg = alogg[wpos[wlogg]-2:wpos[wlogg]]
        endif else begin
            dlogg = reverse(dlogg[wpos[wlogg]-1:wpos[wlogg]+1])
            tlogg = alogg[wpos[wlogg]-1:wpos[wlogg]+1]
        endelse
        ilogg = 3
    endelse

    dfeh = afeh - feh
    w0 = where(abs(dfeh) lt eps, c0)
    if c0 eq 1 then begin
        dfeh = 1d
        tfeh = afeh[w0]
        ifeh = 1
    endif else begin
        wpos = where(dfeh gt 0)
        junk = min(dfeh[wpos], wfeh)
        if (abs(dfeh[wpos[wfeh]-1]) lt abs(dfeh[wpos[wfeh]]) and wpos[wfeh]-1 gt 0) or (wpos[wfeh] ge nfeh) then begin
            dfeh = reverse(dfeh[wpos[wfeh]-2:wpos[wfeh]])
            tfeh = afeh[wpos[wfeh]-2:wpos[wfeh]]
        endif else begin
            dfeh = reverse(dfeh[wpos[wfeh]-1:wpos[wfeh]+1])
            tfeh = afeh[wpos[wfeh]-1:wpos[wfeh]+1]
        endelse
        ifeh = 3
    endelse

    dalpha = aalpha - alpha
    w0 = where(abs(dalpha) lt eps, c0)
    if c0 eq 1 then begin
        dalpha = 1d
        talpha = aalpha[w0]
        ialpha = 1
    endif else begin
        wpos = where(dalpha gt 0)
        junk = min(dalpha[wpos], walpha)
        if (abs(dalpha[wpos[walpha]-1]) lt abs(dalpha[wpos[walpha]]) and wpos[walpha]-1 gt 0) or (wpos[walpha] ge nalpha) then begin
            dalpha = reverse(dalpha[wpos[walpha]-2:wpos[walpha]])
            talpha = aalpha[wpos[walpha]-2:wpos[walpha]]
        endif else begin
            dalpha = reverse(dalpha[wpos[walpha]-1:wpos[walpha]+1])
            talpha = aalpha[wpos[walpha]-1:wpos[walpha]+1]
        endelse
        ialpha = 3
    endelse

    mooglambda = read_moog_lambda7(n=nmoog, fullres=fullres)
    moogspec = dblarr(nmoog)
    for i=0,iteff-1 do begin
        for j=0,ilogg-1 do begin
            for m=0,ifeh-1 do begin
                for n=0,ialpha-1 do begin
                    i1 = (i+1) mod 3
                    i2 = (i+2) mod 3
                    j1 = (j+1) mod 3
                    j2 = (j+2) mod 3
                    m1 = (m+1) mod 3
                    m2 = (m+2) mod 3
                    n1 = (n+1) mod 3
                    n2 = (n+2) mod 3
                    dt = iteff eq 1 ? 1d : dteff[i1]*dteff[i2] / ((tteff[i]-tteff[i1])*(tteff[i]-tteff[i2]))
                    dg = ilogg eq 1 ? 1d : dlogg[j1]*dlogg[j2] / ((tlogg[j]-tlogg[j1])*(tlogg[j]-tlogg[j2]))
                    df = ifeh eq 1 ? 1d : dfeh[m1]*dfeh[m2] / ((tfeh[m]-tfeh[m1])*(tfeh[m]-tfeh[m2]))
                    da = ialpha eq 1 ? 1d : dalpha[n1]*dalpha[n2] / ((talpha[n]-talpha[n1])*(talpha[n]-talpha[n2]))
                    moogspeci = read_moog_bin7(tteff[i], tlogg[j], tfeh[m], talpha[n], n=nmoog, status=status, fullres=fullres)
                    moogspec += dt*dg*df*da*moogspeci
                    wnan = where(finite(moogspeci) eq 0, cnan)
                    if status ne 1 or cnan ne 0 then begin
                        message, 'Synthetic spectrum not found.'
                    endif
                endfor
            endfor
        endfor
    endfor
    return, moogspec
end
