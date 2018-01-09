;function npar, array
;    return, round((double(array[1])-double(array[0]))/double(array[2])) + 1
;end
;
;
;function makepargrid, array, npar
;    return, dindgen(npar)*array[2] + array[0]
;end


function solution, x1, x2, f1, f2, d1, d2
    dx = x2 - x1
    px = x2 + x1
    df = f2 - f1
    pd = d2 + d1

    a = f2*x1*x1*(x1 - 3d * x2) + x2*(x1*dx*(d2*x1+d1*x2) + f1*x2*(3d * x1 - x2))
    b = x2*(6d * x1*df - d1*(px+x1)*dx) - d2*x1*dx*(px+x2)
    c = dx*(pd*px + d2*x1 + d1*x2) - 3d * df*px
    d = 2*df - pd*dx

    return, [[a], [b], [c], [d]] / (x1-x2)^3d
end


function interp_moog7_cubic, teff, logg, feh, alpha
    common hashtable, ht, nmoog

    ;teff_grid1 = [3500., 5501., 100.]
    ;teff_grid2 = logg lt 0.5 ? [5600., 6801., 200.] : [5600., 8001., 200.]
    ;logg_grid = [0.0, 5.01, 0.5]
    ;feh_grid = [-5.0, 0.01, 0.1]
    ;alpha_grid = [-0.8, 1.21, 0.1]
    ;
    ;nteff1 = npar(teff_grid1)
    ;nteff2 = npar(teff_grid2)
    ;nteff = nteff1 + nteff2
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
    if ateff[wteff] eq teff then iteff = wteff else begin
        if wteff eq nteff-1 then wteff -= 1
        iteff = [(wteff-1) > 0, wteff, wteff+1, (wteff+2) < (nteff-1)]
    endelse
    nt = n_elements(iteff)

    wlogg = value_locate(alogg, logg)
    if alogg[wlogg] eq logg then ilogg = wlogg else begin
        if wlogg eq nlogg-1 then wlogg -= 1
        ilogg = [(wlogg-1) > 0, wlogg, wlogg+1, (wlogg+2) < (nlogg-1)]
    endelse
    ng = n_elements(ilogg)

    wfeh = value_locate(afeh, feh)
    if afeh[wfeh] eq feh then ifeh = wfeh else begin
        if wfeh eq nfeh-1 then wfeh -= 1
        ifeh = [(wfeh-1) > 0, wfeh, wfeh+1, (wfeh+2) < (nfeh-1)]
    endelse
    nf = n_elements(ifeh)

    walpha = value_locate(aalpha, alpha)
    if aalpha[walpha] eq alpha then ialpha = walpha else begin
        if walpha eq nalpha-1 then walpha -= 1
        ialpha = [(walpha-1) > 0, walpha, walpha+1, (walpha+2) < (nalpha-1)]
    endelse
    na = n_elements(ialpha)

    if nt gt 1 then steff = dblarr(nmoog,2)
    if ng gt 1 then slogg = dblarr(nmoog,2)
    if nf gt 1 then sfeh = dblarr(nmoog,2)
    if na gt 1 then salpha = dblarr(nmoog,2)

    for i=0,nt-1 do begin
        for j=0,ng-1 do begin
            if iteff[i] ge 28 and ilogg[j] eq 0 then iilogg = 1 else iilogg = ilogg[j]
            for m=0,nf-1 do begin
                for n=0,na-1 do begin
                    keyname = string(iteff[i], iilogg, ifeh[m], ialpha[n], format='(4(I02))')
                    if na eq 1 then begin
                        a = ht->get(keyname)
                        if a[0] eq -1 then begin
                            a = read_moog_bin7(ateff[iteff[i]], alogg[iilogg], afeh[ifeh[m]], aalpha[ialpha[n]], n=nmoog, status=status, fullres=fullres)
                            ht->add, keyname, a
                        endif
                        a = [[a], [a*0d], [a*0d], [a*0d]]
                    endif else begin
                        case n of
                            0: begin
                                dalpha1 = ht->get(keyname)
                                if dalpha1[0] eq -1 then begin
                                    dalpha1 = read_moog_bin7(ateff[iteff[i]], alogg[iilogg], afeh[ifeh[m]], aalpha[ialpha[n]], n=nmoog, status=status, fullres=fullres)
                                    ht->add, keyname, dalpha1
                                endif
                            end
                            1: begin
                                salphatemp = ht->get(keyname)
                                if salphatemp[0] eq -1 then begin
                                    salphatemp = read_moog_bin7(ateff[iteff[i]], alogg[iilogg], afeh[ifeh[m]], aalpha[ialpha[n]], n=nmoog, status=status, fullres=fullres)
                                    ht->add, keyname, salphatemp
                                endif
                                salpha[*,0] = salphatemp
                            end
                            2: begin
                                salphatemp = ht->get(keyname)
                                if salphatemp[0] eq -1 then begin
                                    salphatemp = read_moog_bin7(ateff[iteff[i]], alogg[iilogg], afeh[ifeh[m]], aalpha[ialpha[n]], n=nmoog, status=status, fullres=fullres)
                                    ht->add, keyname, salphatemp
                                endif
                                salpha[*,1] = salphatemp
                                dalpha1 = (salpha[*,1] - dalpha1)/(aalpha[ialpha[2]]-aalpha[ialpha[0]])
                            end
                            3: begin
                                dalpha2 = ht->get(keyname)
                                if dalpha2[0] eq -1 then begin
                                    dalpha2 = read_moog_bin7(ateff[iteff[i]], alogg[iilogg], afeh[ifeh[m]], aalpha[ialpha[n]], n=nmoog, status=status, fullres=fullres)
                                    ht->add, keyname, dalpha2
                                endif
                                dalpha2 = (dalpha2 - salpha[*,0])/(aalpha[ialpha[3]]-aalpha[ialpha[1]])
                            end
                        endcase
                    endelse
                endfor
                if na gt 1 then a = solution(aalpha[ialpha[1]], aalpha[ialpha[2]], salpha[*,0], salpha[*,1], dalpha1, dalpha2)
                if nf eq 1 then begin
                    b = a[*,0] + a[*,1]*alpha + a[*,2]*alpha^2d + a[*,3]*alpha^3d
                    b = [[b], [b*0d], [b*0d], [b*0d]]
                endif else begin
                    case m of
                        0: dfeh1 = a[*,0] + a[*,1]*alpha + a[*,2]*alpha^2d + a[*,3]*alpha^3d
                        1: sfeh[*,m-1] = a[*,0] + a[*,1]*alpha + a[*,2]*alpha^2d + a[*,3]*alpha^3d
                        2: begin
                            sfeh[*,m-1] = a[*,0] + a[*,1]*alpha + a[*,2]*alpha^2d + a[*,3]*alpha^3d
                            dfeh1 = (sfeh[*,m-1] - dfeh1)/(afeh[ifeh[2]]-afeh[ifeh[0]])
                        end
                        3: dfeh2 = (a[*,0] + a[*,1]*alpha + a[*,2]*alpha^2d + a[*,3]*alpha^3d - sfeh[*,0])/(afeh[ifeh[3]]-afeh[ifeh[1]])
                    endcase
                endelse
            endfor
            if nf gt 1 then b = solution(afeh[ifeh[1]], afeh[ifeh[2]], sfeh[*,0], sfeh[*,1], dfeh1, dfeh2)
            if ng eq 1 then begin
                c = b[*,0] + b[*,1]*feh + b[*,2]*feh^2d + b[*,3]*feh^3d
                c = [[c], [c*0d], [c*0d], [c*0d]]
            endif else begin
                case j of
                    0: dlogg1 = b[*,0] + b[*,1]*feh + b[*,2]*feh^2d + b[*,3]*feh^3d
                    1: slogg[*,j-1] = b[*,0] + b[*,1]*feh + b[*,2]*feh^2d + b[*,3]*feh^3d
                    2: begin
                        slogg[*,j-1] = b[*,0] + b[*,1]*feh + b[*,2]*feh^2d + b[*,3]*feh^3d
                        dlogg1 = (slogg[*,j-1] - dlogg1)/(alogg[ilogg[2]]-alogg[ilogg[0]])
                    end
                    3: dlogg2 = (b[*,0] + b[*,1]*feh + b[*,2]*feh^2d + b[*,3]*feh^3d - slogg[*,0])/(alogg[ilogg[3]]-alogg[ilogg[1]])
                endcase
            endelse
        endfor
        if ng gt 1 then c = solution(alogg[ilogg[1]], alogg[ilogg[2]], slogg[*,0], slogg[*,1], dlogg1, dlogg2)
        if nt eq 1 then begin
            d = c[*,0] + c[*,1]*logg + c[*,2]*logg^2d + c[*,3]*logg^3d
            d = [[d], [d*0d], [d*0d], [d*0d]]
        endif else begin
            case i of
                0: dteff1 = c[*,0] + c[*,1]*logg + c[*,2]*logg^2d + c[*,3]*logg^3d
                1: steff[*,i-1] = c[*,0] + c[*,1]*logg + c[*,2]*logg^2d + c[*,3]*logg^3d
                2: begin
                    steff[*,i-1] = c[*,0] + c[*,1]*logg + c[*,2]*logg^2d + c[*,3]*logg^3d
                    dteff1 = (steff[*,i-1] - dteff1)/(ateff[iteff[2]]-ateff[iteff[0]])
                end
                3: dteff2 = (c[*,0] + c[*,1]*logg + c[*,2]*logg^2d + c[*,3]*logg^3d - steff[*,0])/(ateff[iteff[3]]-ateff[iteff[1]])
            endcase
        endelse
    endfor
    if nt gt 1 then d = solution(ateff[iteff[1]], ateff[iteff[2]], steff[*,0], steff[*,1], dteff1, dteff2)
    moogspec = d[*,0] + d[*,1]*teff + d[*,2]*teff^2d + d[*,3]*teff^3d

    htcount = ht->count()
    if htcount ge 3967 then begin
        htkeys = ht->keys()
        for i=0,2047 do ht->remove, htkeys[i]
    endif

    return, moogspec
end



