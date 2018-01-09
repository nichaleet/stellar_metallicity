function solution, x1, x2, f1, f2, d1, d2
    a = f2*x1^2*(x1-3*x2) + x2*(d2*x1^2*(x2-x1) + x2*(f1*(3*x1-x2) + d1*x1*(x2-x1)))
    b = d2*x1*(x1^2 + x1*x2 - 2*x2^2) - x2*(6*x1*(f1-f2) + d1*(x2^2 + x1*x2 - 2*x1^2))
    c = 3*(x1+x2)*(f1-f2) - (x1-x2)*(d1*x1 + 2*d2*x1 + 2*d1*x2 + d2*x2)
    d = 2*(f2-f1) + (d1+d2)*(x1-x2)
    return, [[a], [b], [c], [d]] / (x1-x2)^3
end


function npar, array
    return, round((double(array[1])-double(array[0]))/double(array[2])) + 1
end


function makepargrid, array, npar
    return, dindgen(npar)*array[2] + array[0]
end


function interp_moog7, teff, logg, feh, alpha, mooglambda=mooglambda, fullres=fullres
    teff_grid1 = [3500., 5501., 100.]
    teff_grid2 = logg lt 0.5 ? [5600., 6801., 200.] : [5600., 8001., 200.]
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

    wteff = value_locate(ateff, teff)
    if wteff eq nteff-1 then wteff -= 1
    iteff = [(wteff-1) > 0, wteff, wteff+1, (wteff+2) < (nteff-1)]

    wlogg = value_locate(alogg, logg)
    if wlogg eq nlogg-1 then wlogg -= 1
    ilogg = [(wlogg-1) > 0, wlogg, wlogg+1, (wlogg+2) < (nlogg-1)]

    wfeh = value_locate(afeh, feh)
    if wfeh eq nfeh-1 then wfeh -= 1
    ifeh = [(wfeh-1) > 0, wfeh, wfeh+1, (wfeh+2) < (nfeh-1)]

    walpha = value_locate(aalpha, alpha)
    if walpha eq nalpha-1 then walpha -= 1
    ialpha = [(walpha-1) > 0, walpha, walpha+1, (walpha+2) < (nalpha-1)]

    mooglambda = read_moog_lambda7(n=nmoog, fullres=fullres)
    steff = dblarr(nmoog,2)
    slogg = dblarr(nmoog,2)
    sfeh = dblarr(nmoog,2)
    salpha = dblarr(nmoog,2)
    for i=0,3 do begin
        for j=0,3 do begin
            for m=0,3 do begin
                for n=0,3 do begin
                    case n of
                        0: dalpha1 = read_moog_bin7(ateff[iteff[i]], alogg[ilogg[j]], afeh[ifeh[m]], aalpha[ialpha[n]], n=nmoog, status=status, fullres=fullres)
                        1: salpha[*,n-1] = read_moog_bin7(ateff[iteff[i]], alogg[ilogg[j]], afeh[ifeh[m]], aalpha[ialpha[n]], n=nmoog, status=status, fullres=fullres)
                        2: begin
                            salpha[*,n-1] = read_moog_bin7(ateff[iteff[i]], alogg[ilogg[j]], afeh[ifeh[m]], aalpha[ialpha[n]], n=nmoog, status=status, fullres=fullres)
                            dalpha1 = (salpha[*,n-1] - dalpha1)/(aalpha[ialpha[2]]-aalpha[ialpha[0]])
                        end
                        3: dalpha2 = (read_moog_bin7(ateff[iteff[i]], alogg[ilogg[j]], afeh[ifeh[m]], aalpha[ialpha[n]], n=nmoog, status=status, fullres=fullres) - salpha[*,0])/(aalpha[ialpha[3]]-aalpha[ialpha[1]])
                    endcase
                    ;wnan = where(finite(salpha[*,n]) eq 0, cnan)
                    ;if status ne 1 or cnan ne 0 then begin
                    ;    message, 'Synthetic spectrum not found.'
                    ;endif
                endfor
                c = solution(aalpha[ialpha[1]], aalpha[ialpha[2]], salpha[*,0], salpha[*,1], dalpha1, dalpha2)
                case m of
                    0: dfeh1 = c[*,0] + c[*,1]*alpha + c[*,2]*alpha^2d + c[*,3]*alpha^3d
                    1: sfeh[*,m-1] = c[*,0] + c[*,1]*alpha + c[*,2]*alpha^2d + c[*,3]*alpha^3d
                    2: begin
                        sfeh[*,m-1] = c[*,0] + c[*,1]*alpha + c[*,2]*alpha^2d + c[*,3]*alpha^3d
                        dfeh1 = (sfeh[*,m-1] - dfeh1)/(afeh[ifeh[2]]-afeh[ifeh[0]])
                    end
                    3: dfeh2 = (c[*,0] + c[*,1]*alpha + c[*,2]*alpha^2d + c[*,3]*alpha^3d - sfeh[*,0])/(afeh[ifeh[3]]-afeh[ifeh[1]])
                endcase
            endfor
            c = solution(afeh[ifeh[1]], afeh[ifeh[2]], sfeh[*,0], sfeh[*,1], dfeh1, dfeh2)
            case j of
                0: dlogg1 = c[*,0] + c[*,1]*feh + c[*,2]*feh^2d + c[*,3]*feh^3d
                1: slogg[*,j-1] = c[*,0] + c[*,1]*feh + c[*,2]*feh^2d + c[*,3]*feh^3d
                2: begin
                    slogg[*,j-1] = c[*,0] + c[*,1]*feh + c[*,2]*feh^2d + c[*,3]*feh^3d
                    dlogg1 = (slogg[*,j-1] - dlogg1)/(alogg[ilogg[2]]-alogg[ilogg[0]])
                end
                3: dlogg2 = (c[*,0] + c[*,1]*feh + c[*,2]*feh^2d + c[*,3]*feh^3d - slogg[*,0])/(alogg[ilogg[3]]-alogg[ilogg[1]])
            endcase
        endfor
        c = solution(alogg[ilogg[1]], alogg[ilogg[2]], slogg[*,0], slogg[*,1], dlogg1, dlogg2)
        case i of
            0: dteff1 = c[*,0] + c[*,1]*logg + c[*,2]*logg^2d + c[*,3]*logg^3d
            1: steff[*,i-1] = c[*,0] + c[*,1]*logg + c[*,2]*logg^2d + c[*,3]*logg^3d
            2: begin
                steff[*,i-1] = c[*,0] + c[*,1]*logg + c[*,2]*logg^2d + c[*,3]*logg^3d
                dteff1 = (steff[*,i-1] - dteff1)/(ateff[iteff[2]]-ateff[iteff[0]])
            end
            3: dteff2 = (c[*,0] + c[*,1]*logg + c[*,2]*logg^2d + c[*,3]*logg^3d - steff[*,0])/(ateff[iteff[3]]-ateff[iteff[1]])
        endcase
    endfor
    c = solution(ateff[iteff[1]], ateff[iteff[2]], steff[*,0], steff[*,1], dteff1, dteff2)
    moogspec = c[*,0] + c[*,1]*teff + c[*,2]*teff^2d + c[*,3]*teff^3d
    
    return, moogspec
end



