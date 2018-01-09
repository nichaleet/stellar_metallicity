function interp_moog7, teff, logg, feh, alpha
    common hashtable, ht, nmoog

    ;teff_grid1 = [3500., 5501., 100.]
    ;teff_grid2 = [5600., 8001., 200.]
    ;logg_grid = [0.0, 5.01, 0.5]
    ;feh_grid = [-5.0, 0.01, 0.1]
    ;alpha_grid = [-0.8, 1.21, 0.1]
    ;
    ;nteff1 = npar(teff_grid1)
    ;nteff2 = npar(teff_grid2)
    ;nlogg = npar(logg_grid)
    ;nfeh = npar(feh_grid)
    ;nalpha = npar(alpha_grid)
    ;
    ;ateff = [makepargrid(teff_grid1, nteff1), makepargrid(teff_grid2, nteff2)]
    ;alogg = makepargrid(logg_grid, nlogg)
    ;afeh = makepargrid(feh_grid, nfeh)
    ;aalpha = makepargrid(alpha_grid, nalpha)

    ateff = logg lt 0.5 ? double([3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500, 5600, 5800, 6000, 6200, 6400, 6600, 6800]) : double([3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500, 5600, 5800, 6000, 6200, 6400, 6600, 6800, 7000, 7200, 7400, 7600, 7800, 8000])
    alogg = double([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
    afeh = double([-5.0, -4.9, -4.8, -4.7, -4.6, -4.5, -4.4, -4.3, -4.2, -4.1, -4.0, -3.9, -3.8, -3.7, -3.6, -3.5, -3.4, -3.3, -3.2, -3.1, -3.0, -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0])
    aalpha = double([-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2])

    nteff = n_elements(ateff)
    nlogg = n_elements(alogg)
    nfeh = n_elements(afeh)
    nalpha = n_elements(aalpha)

    wteff = value_locate(ateff, teff)
    if ateff[wteff] eq teff then begin
        iteff = wteff
        dteff = 1d
        nt = 1
    endif else begin
        if wteff eq nteff-1 then wteff -= 1
        iteff = [wteff, wteff+1]
        dteff = abs(reverse(ateff[iteff])-teff)
        nt = 2
    endelse

    wlogg = value_locate(alogg, logg)
    if alogg[wlogg] eq logg then begin
        ilogg = wlogg
        dlogg = 1d
        ng = 1
    endif else begin
        if wlogg eq nlogg-1 then wlogg -= 1
        ilogg = [wlogg, wlogg+1]
        dlogg = abs(reverse(alogg[ilogg])-logg)
        ng = 2
    endelse

    wfeh = value_locate(afeh, feh)
    if afeh[wfeh] eq feh then begin
        ifeh = wfeh
        dfeh = 1d
        nf = 1
    endif else begin
        if wfeh eq nfeh-1 then wfeh -= 1
        ifeh = [wfeh, wfeh+1]
        dfeh = abs(reverse(afeh[ifeh])-feh)
        nf = 2
    endelse

    walpha = value_locate(aalpha, alpha)
    if aalpha[walpha] eq alpha then begin
        ialpha = walpha
        dalpha = 1d
        na = 1
    endif else begin
        if walpha eq nalpha-1 then walpha -= 1
        ialpha = [walpha, walpha+1]
        dalpha = abs(reverse(aalpha[ialpha])-alpha)
        na = 2
    endelse

    moogspec = dblarr(nmoog)
    for i=0,nt-1 do begin
        for j=0,ng-1 do begin
            for m=0,nf-1 do begin
                for n=0,na-1 do begin
                    keyname = string(iteff[i], ilogg[j], ifeh[m], ialpha[n], format='(4(I02))')
                    moogspeci = ht->get(keyname)
                    if moogspeci[0] eq -1 then begin
                        moogspeci = read_moog_bin7(ateff[iteff[i]], alogg[ilogg[j]], afeh[ifeh[m]], aalpha[ialpha[n]], n=nmoog)
                        ht->add, keyname, moogspeci
                    endif
                    moogspec += dteff[i]*dlogg[j]*dfeh[m]*dalpha[n]*moogspeci
                endfor
            endfor
        endfor
    endfor
    moogspec /= (nt eq 2 ? (ateff[iteff[1]]-ateff[iteff[0]]) : 1d) $
                *(ng eq 2 ? (alogg[ilogg[1]]-alogg[ilogg[0]]) : 1d) $
                *(nf eq 2 ? (afeh[ifeh[1]]-afeh[ifeh[0]]) : 1d) $
                *(na eq 2 ? (aalpha[ialpha[1]]-aalpha[ialpha[0]]) : 1d)

    return, moogspec
end
