function npar, array
    return, round((double(array[1])-double(array[0]))/double(array[2])) + 1
end


function makepargrid, array, npar
    return, dindgen(npar)*array[2] + array[0]
end


function interp_moog, teff, logg, feh, alpha, vt, mooglambda=mooglambda, fullres=fullres
    teff_grid1 = [4000., 5501., 100.]
    teff_grid2 = [5600., 8001., 200.]
    logg_grid = [0.0, 3.51, 0.5]
    feh_grid = [-4.0, 0.01, 0.1]
    alpha_grid = [-0.6, 0.61, 0.1]
    vt_grid = [0.0, 3.01, 0.5]

    nteff1 = npar(teff_grid1)
    nteff2 = npar(teff_grid2)
    nlogg = npar(logg_grid)
    nfeh = npar(feh_grid)
    nalpha = npar(alpha_grid)
    nvt = npar(vt_grid)

    ateff = [makepargrid(teff_grid1, nteff1), makepargrid(teff_grid2, nteff2)]
    alogg = makepargrid(logg_grid, nlogg)
    afeh = makepargrid(feh_grid, nfeh)
    aalpha = makepargrid(alpha_grid, nalpha)
    avt = makepargrid(vt_grid, nvt)

    dteff = ateff - teff
    w0 = where(abs(dteff) lt 1d-5, c0)
    if c0 eq 1 then begin
        dteff = 1d
        tteff = ateff[w0]
        iteff = 1
    endif else begin
        wpos = where(dteff gt 0)
        junk = min(dteff[wpos], wteff)
        dteff = reverse(abs(dteff[wpos[wteff]-1:wpos[wteff]]))
        tteff = ateff[wpos[wteff]-1:wpos[wteff]]
        iteff = 2
    endelse

    dlogg = alogg - logg
    w0 = where(abs(dlogg) lt 1d-5, c0)
    if c0 eq 1 then begin
        dlogg = 1d
        tlogg = alogg[w0]
        ilogg = 1
    endif else begin
        wpos = where(dlogg gt 0)
        junk = min(dlogg[wpos], wlogg)
        dlogg = reverse(abs(dlogg[wpos[wlogg]-1:wpos[wlogg]]))
        tlogg = alogg[wpos[wlogg]-1:wpos[wlogg]]
        ilogg = 2
    endelse

    dfeh = afeh - feh
    w0 = where(abs(dfeh) lt 1d-5, c0)
    if c0 eq 1 then begin
        dfeh = 1d
        tfeh = afeh[w0]
        ifeh = 1
    endif else begin
        wpos = where(dfeh gt 0)
        junk = min(dfeh[wpos], wfeh)
        dfeh = reverse(abs(dfeh[wpos[wfeh]-1:wpos[wfeh]]))
        tfeh = afeh[wpos[wfeh]-1:wpos[wfeh]]
        ifeh = 2
    endelse

    dalpha = aalpha - alpha
    w0 = where(abs(dalpha) lt 1d-5, c0)
    if c0 eq 1 then begin
        dalpha = 1d
        talpha = aalpha[w0]
        ialpha = 1
    endif else begin
        wpos = where(dalpha gt 0)
        junk = min(dalpha[wpos], walpha)
        dalpha = reverse(abs(dalpha[wpos[walpha]-1:wpos[walpha]]))
        talpha = aalpha[wpos[walpha]-1:wpos[walpha]]
        ialpha = 2
    endelse

    ivt = 2

    mooglambda = read_moog_lambda(n=nmoog, fullres=fullres)
    moogspec = dblarr(nmoog)
    moogspectot = dblarr(nmoog)
    replace = bytarr(iteff, ilogg, ifeh, ialpha, ivt)
    for i=0,iteff-1 do begin
        for j=0,ilogg-1 do begin
            for m=0,ifeh-1 do begin
                for n=0,ialpha-1 do begin
                    vt = 2.7 - 0.509*tlogg[j]
                    dvt = avt - vt
                    wpos = where(dvt gt 0)
                    junk = min(dvt[wpos], wvt)
                    dvt = reverse(abs(dvt[wpos[wvt]-1:wpos[wvt]]))
                    tvt = avt[wpos[wvt]-1:wpos[wvt]]
                    moogspeci = (dvt[0]*read_moog_bin(tteff[i], tlogg[j], tfeh[m], talpha[n], tvt[0], n=nmoog, status=status, fullres=fullres) + $
                                 dvt[1]*read_moog_bin(tteff[i], tlogg[j], tfeh[m], talpha[n], tvt[1], n=nmoog, status=status, fullres=fullres)) / $
                                (tvt[1]-tvt[0])
                    moogspec += dteff[i]*dlogg[j]*dfeh[m]*dalpha[n]*moogspeci
                    moogspectot += moogspeci
                    wnan = where(finite(moogspeci) eq 0, cnan)
                    if status ne 1 or cnan ne 0 then begin
                        ;if cnan gt 0 then print, tteff[i], tlogg[j], tfeh[m], talpha[n], tvt[p], cnan
                        replace[i,j,m,n,p] = 1
                        continue
                    endif
                endfor
            endfor
        endfor
    endfor
    wr = where(replace eq 1, cr, ncomplement=ngood)
    if cr gt 0 then begin
        message, 'Missing grid point!'
        moogspectot /= ngood
        for i=0,cr-1 do begin
            ar = array_indices(replace, wr[i])
            moogspec += dteff[ar[0]]*dlogg[ar[1]]*dfeh[ar[2]]*dalpha[ar[3]]*dalpha[ar[4]]*moogspectot
        endfor
    endif
    moogspec /= (iteff eq 2 ? (tteff[1]-tteff[0]) : 1d) $
                *(ilogg eq 2 ? (tlogg[1]-tlogg[0]) : 1d) $
                *(ifeh eq 2 ? (tfeh[1]-tfeh[0]) : 1d) $
                *(ialpha eq 2 ? (talpha[1]-talpha[0]) : 1d)

    return, moogspec
end
