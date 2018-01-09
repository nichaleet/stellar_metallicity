function interp_moog_ch, teff, logg, feh, alpha, cfe, fullres=fullres, mds=mds
    common hashtable_ch, ht_ch, nmoog_ch

    ateff = double([3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500, 5600, 5800, 6000, 6200, 6400])
    alogg = double([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0])
    afeh = double([-4.0, -3.9, -3.8, -3.7, -3.6, -3.5, -3.4, -3.3, -3.2, -3.1, -3.0, -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0])
    aalpha = double([-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2])
    acfe = double([-2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2,  0.0,  0.2, 0.4,  0.6,  0.8,  1.0,  1.4, 1.8,  2.2,  2.6,  3.0,  3.5])

    nteff = n_elements(ateff)
    nlogg = n_elements(alogg)
    nfeh = n_elements(afeh)
    nalpha = n_elements(aalpha)
    ncfe = n_elements(acfe)

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

    wcfe = value_locate(acfe, cfe)
    if acfe[wcfe] eq cfe then begin
        icfe = wcfe
        dcfe = 1d
        nc = 1
    endif else begin
        if wcfe eq ncfe-1 then wcfe -= 1
        icfe = [wcfe, wcfe+1]
        dcfe = abs(reverse(acfe[icfe])-cfe)
        nc = 2
    endelse

    moogspec_ch = dblarr(nmoog_ch)
    for i=0,nt-1 do begin
        for j=0,ng-1 do begin
            for m=0,nf-1 do begin
                for n=0,na-1 do begin
                    for p=0,nc-1 do begin
                        keyname = string(iteff[i], ilogg[j], ifeh[m], ialpha[n], icfe[p], format='(5(I02))')
                        moogspeci_ch = ht_ch->get(keyname)
                        if moogspeci_ch[0] eq -1 then begin
                            moogspeci_ch = read_moog_bin_ch(ateff[iteff[i]], alogg[ilogg[j]], afeh[ifeh[m]], aalpha[ialpha[n]], acfe[icfe[p]], fullres=fullres, mds=mds)
                            ht_ch->add, keyname, moogspeci_ch
                        endif
                        moogspec_ch += dteff[i]*dlogg[j]*dfeh[m]*dalpha[n]*dcfe[p]*moogspeci_ch
                    endfor
                endfor
            endfor
        endfor
    endfor
    moogspec_ch /= (nt eq 2 ? (ateff[iteff[1]]-ateff[iteff[0]]) : 1d) $
                *(ng eq 2 ? (alogg[ilogg[1]]-alogg[ilogg[0]]) : 1d) $
                *(nf eq 2 ? (afeh[ifeh[1]]-afeh[ifeh[0]]) : 1d) $
                *(na eq 2 ? (aalpha[ialpha[1]]-aalpha[ialpha[0]]) : 1d) $
                *(nc eq 2 ? (acfe[icfe[1]]-acfe[icfe[0]]) : 1d)

    return, moogspec_ch
end
