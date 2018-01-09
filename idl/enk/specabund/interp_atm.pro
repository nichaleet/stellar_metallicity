pro interp_atm, teff, logg, vt, feh, alphafe, outfile=outfile, tweakel=tweakel, tweakabund=tweakabund
    if ~keyword_set(outfile) then outfile = 'interp.atm'
    
    path = getenv('chome')+'atlas/atlas9/'
    gridteff = [dindgen(21)*100+3500, dindgen(13)*200+5600]
    gridlogg = dindgen(11)*0.5
    gridvt = [0d, 1d, 2d, 4d]
    gridfeh = dindgen(11)*0.5-5.0
    gridalphafe = dindgen(21)*0.1-0.8

    if teff gt max(gridteff) or teff lt min(gridteff) then message, 'teff out of range.'
    if logg gt max(gridlogg) or logg lt min(gridlogg) then message, 'logg out of range.'
    if feh gt max(gridfeh) or feh lt min(gridfeh) then message, 'feh out of range.'
    if alphafe gt max(gridalphafe) or alphafe lt min(gridalphafe) then message, 'alphafe out of range.', /info
    alphafe >= min(gridalphafe)
    alphafe <= max(gridalphafe)
    if vt lt 0 then message, 'vt out of range.'

    dteff = gridteff - teff
    w = where(dteff lt 0)
    junk = max(dteff[w], m)
    m = w[m]
    ateff = gridteff[[m, m+1]]
    dteff = reverse(dteff[[m, m+1]])
    dteff = abs(dteff/(max(dteff)-min(dteff)))

    if logg eq 0.0 then begin
        alogg = [0.0, 0.5]
        dlogg = [1.0, 0.0]
    endif else begin
        dlogg = gridlogg - logg
        w = where(dlogg lt 0)
        junk = max(dlogg[w], m)
        m = w[m]
        alogg = gridlogg[[m, m+1]]
        dlogg = reverse(dlogg[[m, m+1]])
        dlogg = abs(dlogg/(max(dlogg)-min(dlogg)))
    endelse

    dfeh = gridfeh - feh
    w = where(dfeh lt 0)
    junk = max(dfeh[w], m)
    m = w[m]
    afeh = gridfeh[[m, m+1]]
    dfeh = reverse(dfeh[[m, m+1]])
    dfeh = abs(dfeh/(max(dfeh)-min(dfeh)))

    dalphafe = gridalphafe - alphafe
    w = where(dalphafe lt 0)
    junk = max(dalphafe[w], m)
    m = w[m]
    aalphafe = gridalphafe[[m, m+1]]
    dalphafe = reverse(dalphafe[[m, m+1]])
    dalphafe = abs(dalphafe/(max(dalphafe)-min(dalphafe)))

    c1 = dblarr(72)
    c2 = dblarr(72)
    c3 = dblarr(72)
    c4 = dblarr(72)
    c5 = dblarr(72)
    c6 = dblarr(72)
    c7 = dblarr(72)

    for i=0,1 do begin
        for j=0,1 do begin
            for k=0,1 do begin
                for l=0,1 do begin
                    vtcalc = 2.1427009 - 0.23740876*alogg[l]
                    dvt = gridvt - vtcalc
                    w = where(dvt lt 0)
                    junk = max(dvt[w], m)
                    avt = gridvt[[m, m+1]]
                    dvt = reverse(dvt[[m, m+1]])
                    dvt = abs(dvt/(max(dvt)-min(dvt)))
                    case 1 of
                        vt lt avt[0]: begin
                            avt = [avt[0], avt[0]]
                            dvt = [1d, 0d]
                        end
                        vt gt avt[1]: begin
                            avt = [avt[1], avt[1]]
                            dvt = [1d, 0d]
                        end
                        vt ge avt[0] and vt le avt[1]:
                    endcase
                    for m=0,1 do begin
                        factor = dfeh[i]*dalphafe[j]*dteff[k]*dlogg[l]*dvt[m]
                        if factor eq 0.0 then continue

                        fsign = afeh[i] lt 0.0 ? 'm' : 'p'
                        sfeh = string(fsign, round(abs(afeh[i]*10)), format='(A1,I02)')
                        asign = aalphafe[j] lt 0.0 ? 'm' : 'p'
                        salphafe = 'a' + asign + strtrim(round(abs(aalphafe[j]*10)), 2)
                        if aalphafe[j] eq 0.0 then salphafe = ''
                        steff = string('t', ateff[k], format='(A1,I4)')
                        slogg = string('g', round(alogg[l]*10), format='(A1,I02)')
                        svt = string('k', avt[m], format='(A1,I1)')
                        atmname = path+sfeh+salphafe+'/a'+sfeh+salphafe+steff+slogg+svt+'odfnew.dat'

                        readcol, atmname, c1_in, c2_in, c3_in, c4_in, c5_in, c6_in, c7_in, format='D,D,D,D,D,D,D', numline=72, skipline=23, /silent
                        c1 += factor*c1_in
                        c2 += factor*c2_in
                        c3 += factor*c3_in
                        c4 += factor*c4_in
                        c5 += factor*c5_in
                        c6 += factor*c6_in
                        c7 += factor*c7_in
                    endfor
                endfor
            endfor
        endfor
    endfor

    solar = [12.00,10.99, 3.31, 1.42, 2.88, 8.56, 8.05, 8.93, 4.56, 8.09, $
             6.33, 7.58, 6.47, 7.55, 5.45, 7.21, 5.5 , 6.56, 5.12, 6.36, $
             3.10, 4.99, 4.00, 5.67, 5.39, 7.52, 4.92, 6.25, 4.21, 4.60, $
             2.88, 3.41, 2.37, 3.35, 2.63, 3.23, 2.60, 2.90, 2.24, 2.60, $
             1.42, 1.92, 0.00, 1.84, 1.12, 1.69, 1.24, 1.86, 0.82, 2.0 , $
             1.04, 2.24, 1.51, 2.23, 1.12, 2.13, 1.22, 1.55, 0.71, 1.50, $
             0.00, 1.00, 0.51, 1.12, 0.33, 1.1 , 0.50, 0.93, 0.13, 1.08, $
             0.12, 0.88, 0.13, 0.68, 0.27, 1.45, 1.35, 1.8 , 0.83, 1.09, $
             0.82, 1.85, 0.71, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.12, $
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] - 12.0
    w = [8, 10, 12, 14, 16, 18, 20, 22]-1
    out = solar
    out[w] += alphafe
    sum = total(10.^(out[2:n_elements(out)-1]+feh))
    hfrac = 0.92040 / (0.92040 + 0.07834)
    hefrac = 1.0 - hfrac
    h = (1.0 - sum)*hfrac
    he = (1.0 - sum)*hefrac
    out[0] = h
    out[1] = he

    openw, lun, outfile, /get_lun
    ;printf, lun, teff, logg, format='("TEFF",3X,D5.0,2X,"GRAVITY",1X,D7.5,1X,"LTE")'
    ;printf, lun, feh, vt, format='("TITLE  [",D4.1,"] VTURB=",D4.2,2X,"L/H=1.25 NOVER NEW ODF")'
    ;printf, lun, ' OPACITY IFOP 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0'
    ;printf, lun, ' CONVECTION ON   1.25 TURBULENCE OFF  0.00  0.00  0.00  0.00'
    ;printf, lun, 10^(feh[f]), out[0], out[1], format='("ABUNDANCE SCALE ",D9.5," ABUNDANCE CHANGE 1 ",D7.5," 2 ",D7.5)'
    ;printf, lun, out[2:7], format='(" ABUNDANCE CHANGE  3 ",D6.2,"  4 ",D6.2,"  5 ",D6.2,"  6 ",D6.2,"  7 ",D6.2,"  8 ",D6.2)'
    ;printf, lun, out[8:13], format='(" ABUNDANCE CHANGE  9 ",D6.2," 10 ",D6.2," 11 ",D6.2," 12 ",D6.2," 13 ",D6.2," 14 ",D6.2)'
    ;printf, lun, out[14:19], format='(" ABUNDANCE CHANGE 15 ",D6.2," 16 ",D6.2," 17 ",D6.2," 18 ",D6.2," 19 ",D6.2," 20 ",D6.2)'
    ;printf, lun, out[20:25], format='(" ABUNDANCE CHANGE 21 ",D6.2," 22 ",D6.2," 23 ",D6.2," 24 ",D6.2," 25 ",D6.2," 26 ",D6.2)'
    ;printf, lun, out[26:31], format='(" ABUNDANCE CHANGE 27 ",D6.2," 28 ",D6.2," 29 ",D6.2," 30 ",D6.2," 31 ",D6.2," 32 ",D6.2)'
    ;printf, lun, out[32:37], format='(" ABUNDANCE CHANGE 33 ",D6.2," 34 ",D6.2," 35 ",D6.2," 36 ",D6.2," 37 ",D6.2," 38 ",D6.2)'
    ;printf, lun, out[38:43], format='(" ABUNDANCE CHANGE 39 ",D6.2," 40 ",D6.2," 41 ",D6.2," 42 ",D6.2," 43 ",D6.2," 44 ",D6.2)'
    ;printf, lun, out[44:49], format='(" ABUNDANCE CHANGE 45 ",D6.2," 46 ",D6.2," 47 ",D6.2," 48 ",D6.2," 49 ",D6.2," 50 ",D6.2)'
    ;printf, lun, out[50:55], format='(" ABUNDANCE CHANGE 51 ",D6.2," 52 ",D6.2," 53 ",D6.2," 54 ",D6.2," 55 ",D6.2," 56 ",D6.2)'
    ;printf, lun, out[56:61], format='(" ABUNDANCE CHANGE 57 ",D6.2," 58 ",D6.2," 59 ",D6.2," 60 ",D6.2," 61 ",D6.2," 62 ",D6.2)'
    ;printf, lun, out[62:67], format='(" ABUNDANCE CHANGE 63 ",D6.2," 64 ",D6.2," 65 ",D6.2," 66 ",D6.2," 67 ",D6.2," 68 ",D6.2)'
    ;printf, lun, out[68:73], format='(" ABUNDANCE CHANGE 69 ",D6.2," 70 ",D6.2," 71 ",D6.2," 72 ",D6.2," 73 ",D6.2," 74 ",D6.2)'
    ;printf, lun, out[74:79], format='(" ABUNDANCE CHANGE 75 ",D6.2," 76 ",D6.2," 77 ",D6.2," 78 ",D6.2," 79 ",D6.2," 80 ",D6.2)'
    ;printf, lun, out[80:85], format='(" ABUNDANCE CHANGE 81 ",D6.2," 82 ",D6.2," 83 ",D6.2," 84 ",D6.2," 85 ",D6.2," 86 ",D6.2)'
    ;printf, lun, out[86:91], format='(" ABUNDANCE CHANGE 87 ",D6.2," 88 ",D6.2," 89 ",D6.2," 90 ",D6.2," 91 ",D6.2," 92 ",D6.2)'
    ;printf, lun, out[92:97], format='(" ABUNDANCE CHANGE 93 ",D6.2," 94 ",D6.2," 95 ",D6.2," 96 ",D6.2," 97 ",D6.2," 98 ",D6.2)'
    ;printf, lun, out[98], format='(" ABUNDANCE CHANGE 99 ",D6.2)'
    ;printf, lun, 'READ DECK6 72 RHOX,T,P,XNE,ABROSS,ACCRAD,VTURB, FLXCNV,VCONV,VELSND'

    printf, lun, 'KURUCZ'
    printf, lun, teff, logg, feh, alphafe, vt, format='(D5.0,"/",D4.2,"/",D+5.2,"/",D+5.2,"/",D4.2)'
    printf, lun, 'ntau=      72'
    for i=0,n_elements(c1)-1 do printf, lun, c1[i], c2[i], c3[i], c4[i], c5[i], c6[i], vt, format='(1X,E15.9,2X,F8.1,5(1X,E10.4))'

    ;printf, lun, 'PRADK 5.9467E-01'
    ;printf, lun, 'BEGIN                    ITERATION  15 COMPLETED'
    
    printf, lun, vt, format='(E13.3)'

    ntweak = 0
    tweakels = [0]
    tweakabunds = [0.0]
    if alphafe ne 0.0 then begin
        ntweak += 8
        tweakels = [tweakels, 8, 10, 12, 14, 16, 17, 20, 22]
        tweakabunds = [tweakabunds, replicate(alphafe, 8)]
    endif
    if keyword_set(tweakel) then begin
        ntweak += n_elements(tweakel)
        tweakels = [tweakels, tweakel]
        tweakabunds = [tweakabunds, tweakabund]
    endif
    
    printf, lun, ntweak, feh, format='("NATOMS",4X,I2,2X,D6.2)'
    s = sort(tweakels)
    for i=0,ntweak-1 do begin
        abund = tweakels[s[i+1]] eq 3 ? tweakabunds[s[i+1]] : solar[tweakels[s[i+1]]-1]+feh+tweakabunds[s[i+1]]+12.0
        printf, lun, tweakels[s[i+1]], abund, format='("      ",I2,"    ",D6.2)'
    endfor

    printf, lun, 'NMOL       18'
    printf, lun, '101.0   106.0   107.0   108.0   606.0   607.0   608.0   707.0'
    printf, lun, '708.0   808.0 10108.0 60808.0     6.1     7.1     8.1    22.1'
    printf, lun, ' 23.1   823.0'
    close, lun
    free_lun, lun
end
