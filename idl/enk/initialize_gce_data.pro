pro initialize_gce_data
    common SN, z_II, nel, nom06, kar10, eps_sun, M_SN, M_HN, M_lims, z_lims, perets

    ;readcol, 'stellar_lifetimes.txt', sn_age, sn_mass, format='D,D', comment='#', /silent

    z_II = [0.0, 0.001, 0.004, 0.02]
    M_SN = [13.0, 15.0, 18.0, 20.0, 25.0, 30.0, 40.0]
    M_HN = [20.0, 25.0, 30.0, 40.0]
    Mfinal_SN = dblarr(4,7)
    Mcut_SN = dblarr(4,7)
    Mfinal_HN = dblarr(4,4)
    Mcut_HN = dblarr(4,4)

    nom06sn = {isotope:0, atomic:0, mass:0d, II:dblarr(4,7), HN:dblarr(4,4), Ia:0d}
    nom06sn = replicate(nom06sn, 100)
    n = 0L
    openr, 1, getenv('M31')+'papers/nom06/tab1.txt'  ;SN II yields
    skip_lun, 1, 24, /lines
    isotope = ' '
    while ~eof(1) do begin
        readf, 1, II_z, isotope, II_13, II_15, II_18, II_20, II_25, II_30, II_40, format='(D5,1X,A8,7(1X,E9))'
        w2 = where(z_II eq II_z)
        if strtrim(isotope, 2) eq 'M_final_' then Mfinal_SN[w2,*] = [II_13, II_15, II_18, II_20, II_25, II_30, II_40] else if strtrim(isotope, 2) eq 'M_cut_' then Mcut_SN[w2,*] = [II_13, II_15, II_18, II_20, II_25, II_30, II_40] else begin
            isotope = strtrim(isotope, 2)
            isotope_AZ = isotope_to_AZ(isotope)
            if II_z eq 0.0 then begin
                w1 = n
                n++
            endif else w1 = where(nom06sn.isotope eq isotope_AZ[0] and nom06sn.atomic eq isotope_AZ[1])
            nom06sn[w1].isotope = isotope_AZ[0]
            nom06sn[w1].atomic = isotope_AZ[1]
            nom06sn[w1].II[w2,*] = [II_13, II_15, II_18, II_20, II_25, II_30, II_40]
        endelse
    endwhile
    close, 1
    nom06sn = nom06sn[0:n-1]

    openr, 1, getenv('M31')+'papers/nom06/tab2.txt'  ;HN II yields
    skip_lun, 1, 21, /lines
    isotope = ' '
    while ~eof(1) do begin
        readf, 1, II_z, isotope, II_20, II_25, II_30, II_40, format='(D5,1X,A8,4(1X,E9))'
        w2 = where(z_II eq II_z)
        if strtrim(isotope, 2) eq 'M_final_' then Mfinal_HN[w2,*] = [II_20, II_25, II_30, II_40] else if strtrim(isotope, 2) eq 'M_cut_' then Mcut_HN[w2,*] = [II_20, II_25, II_30, II_40] else begin
            isotope = strtrim(isotope, 2)
            isotope_AZ = isotope_to_AZ(isotope)
            w1 = where(nom06sn.isotope eq isotope_AZ[0] and nom06sn.atomic eq isotope_AZ[1])
            nom06sn[w1].HN[w2,*] = [II_20, II_25, II_30, II_40]
        endelse
    endwhile
    close, 1

    nom06iso = {isotope:0, atomic:0, II:dblarr(4), Ia:0d}
    nom06iso = replicate(nom06iso, 100)
    n = 0L
    openr, 1, getenv('M31')+'papers/nom06/tab3.txt'  ;SN II + HN yields integrated over IMF
    skip_lun, 1, 16, /lines
    while ~eof(1) do begin
        readf, 1, isotope, II_z0, II_z001, II_z004, II_z02, Ia, format='(A6,5(1X,E8))'
        isotope = strtrim(isotope, 2)
        isotope_AZ = isotope_to_AZ(isotope)
        nom06iso[n].isotope = isotope_AZ[0]
        nom06iso[n].atomic = isotope_AZ[1]
        nom06iso[n].II[0] = II_z0
        nom06iso[n].II[1] = II_z001
        nom06iso[n].II[2] = II_z004
        nom06iso[n].II[3] = II_z02
        nom06iso[n].Ia = Ia
        w = where(nom06sn.isotope eq isotope_AZ[0] and nom06sn.atomic eq isotope_AZ[1])
        nom06sn[w].Ia = Ia
        n++
    endwhile
    close, 1
    nom06iso = nom06iso[0:n-1]

    openr, 1, getenv('M31')+'papers/iwa99/tab3.txt'  ;SN Ia yields (Iwamoto et al. 1999)
    skip_lun, 1, 1, /lines
    while ~eof(1) do begin
        readf, 1, isotope, Ia, format='(A6,13X,E8)'
        isotope = strtrim(isotope, 2)
        isotope_AZ = isotope_to_AZ(isotope)
        w = where(nom06sn.isotope eq isotope_AZ[0] and nom06sn.atomic eq isotope_AZ[1])
        nom06sn[w].Ia = Ia
        n++
    endwhile
    close, 1

    atomic_weight = [1.00794, 4.00602, 6.941, 9.012182, 10.811, 12.0107, 14.0067, 15.9994, 18.9984032, 20.1797, 22.98976928, 24.3050, 26.9815386, 28.0355, 30.973762, 32.065, 35.453, 39.948, 39.0983, 40.078, 44.955912, 47.957, 50.9415, 51.9951, 54.9438045, 55.845, 58.933195, 58.6934, 63.546, 65.38, 69.723, 72.64]
    eps_sun = [12.00, 10.99, 3.31, 1.42, 2.88, 8.56, 8.05, 8.93, 4.56, 8.09, 6.33, 7.58, 6.47, 7.55, 5.45, 7.21, 5.5 , 6.56, 5.12, 6.36, 3.10, 4.99, 4.00, 5.67, 5.39, 7.52, 4.92, 6.25, 4.21, 4.60, 2.88, 3.41]
    el = nom06iso[uniq(nom06iso.atomic, sort(nom06iso.atomic))].atomic
    nel = n_elements(el)
    nom06 = {atomic:0, weight_II:dblarr(4,7), II:dblarr(4,7), weight_HN:dblarr(4,4), HN:dblarr(4,4), weight_Ia:0d, Ia:0d}
    nom06 = replicate(nom06, nel)
    nom06i = {atomic:0, weight:dblarr(4), II:dblarr(4), weight_Ia:0d, Ia:0d}
    nom06i = replicate(nom06i, nel)
    for i=0,nel-1 do begin
        w = where(nom06sn.atomic eq el[i])
        nom06[i].atomic = el[i]
        for j=0,3 do begin
            for k=0,6 do begin
                nom06[i].weight_II[j,k] = total(nom06sn[w].isotope*nom06sn[w].II[j,k])/total(nom06sn[w].II[j,k])
                nom06[i].II[j,k] = total(nom06sn[w].II[j,k])
            endfor
        endfor
        for j=0,3 do begin
            for k=0,3 do begin
                nom06[i].weight_HN[j,k] = total(nom06sn[w].isotope*nom06sn[w].HN[j,k])/total(nom06sn[w].HN[j,k])
                nom06[i].HN[j,k] = total(nom06sn[w].HN[j,k])
            endfor
        endfor
        nom06[i].weight_Ia = total(nom06sn[w].isotope*nom06sn[w].Ia)/total(nom06sn[w].Ia)
        nom06[i].Ia = total(nom06sn[w].Ia)

        w = where(nom06iso.atomic eq el[i])
        nom06i[i].atomic = el[i]
        for j=0,3 do begin
            nom06i[i].weight[j] = total(nom06iso[w].isotope*nom06iso[w].II[j])/total(nom06iso[w].II[j])
            nom06i[i].II[j] = total(nom06iso[w].II[j])
        endfor
        nom06i[i].weight_Ia = total(nom06iso[w].isotope*nom06iso[w].Ia)/total(nom06iso[w].Ia)
        nom06i[i].Ia = total(nom06iso[w].Ia)
    endfor

    M_lims = [1.00, 1.25, 1.50, 1.75, 1.90, 2.00, 2.10, 2.25, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50]
    z_lims = [0.0001, 0.004, 0.008, 0.02]
    kar10files = getenv('M31')+'papers/kar10/table_a'+['2', '3', '4', '5', '6']+'.txt'
    nkar10 = n_elements(kar10files)
    kar10 = {atomic:0, masslost:dblarr(4,17)-999}
    kar10 = replicate(kar10, 19)
    kar10.atomic = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 26, 27, 28]
    iso_in = ' '
    for i=0,nkar10-1 do begin
        openr, 1, kar10files[i]
        while ~eof(1) do begin
            readf, 1, m_init, z_in, format='(13X,D5,11X,D6)'
            wm = where(M_lims eq m_init)
            wz = where(z_lims eq z_in)
            skip_lun, 1, 3, /lines
            for j=0,74 do begin
                readf, 1, iso_in, masslost_in, format='(A7,24X,E14)'
                iso_in = strtrim(iso_in, 2)
                iso_ss = strsplit(iso_in, '123456789-*', length=l)
                elname = strmid(iso_in, 0, l[0])
                ;isotope = strlen(iso_in) lt l[0] ? 1 : fix(strmid(iso_in, l[0])
                case elname of
                    'p': atomic = 1
                    'd': begin
                        atomic = 1
                        ;isotope = 2
                    end
                    'he': atomic = 2
                    'li': atomic = 3
                    'be': atomic = 4
                    'b': atomic = 5
                    'c': atomic = 6
                    'n': atomic = 7
                    'o': atomic = 8
                    'f': atomic = 9
                    'ne': atomic = 10
                    'na': atomic = 11
                    'mg': atomic = 12
                    'al': atomic = 13
                    'si': atomic = 14
                    'p': atomic = 15
                    's': atomic = 16
                    'fe': atomic = 26
                    'co': atomic = 27
                    'ni': atomic = 28
                endcase
                wa = where(kar10.atomic eq atomic)
                if kar10[wa].masslost[wz,wm] eq -999 then kar10[wa].masslost[wz,wm] = 0d
                kar10[wa].masslost[wz,wm] += masslost_in
            endfor
        endwhile
        close, 1
    endfor

    for i=0,n_elements(kar10)-1 do begin
        for j=0,n_elements(z_lims)-1 do begin
            w = where(kar10[i].masslost[j,*] eq -999, complement=ww, c)
            if c eq 0 then continue
            if c eq n_elements(M_lims) then begin
                kar10[i].masslost[j,w] = 0d
                continue
            endif
            kar10[i].masslost[j,w] = interpol(kar10[i].masslost[j,ww], M_lims[ww], M_lims[w], /spline)
        endfor
    endfor

    wel = [0, 1, 7, 11, 13, 19, 21, 25]
    nel = n_elements(wel)
    nom06 = nom06[wel]
    nom06i = nom06i[wel]
    eps_sun = eps_sun[wel]

    perets = replicate({dotIa_1:0d, dotIa_2:0d}, 8)
    perets[3:7].dotIa_1 = [4.4d-5, 0.001, 0.07, 0.0059, 0.00025]
    perets[3:7].dotIa_2 = [2.7d-5, 0.00058, 0.052, 0.021, 0.0025]

    wel = [0, 1, 7, 11, 13, 16]
    kar10 = kar10[wel]
end
