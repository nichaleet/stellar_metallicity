pro atlas9_mach, teff, logg, vt, feh, alphafe
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
    c8 = dblarr(72)
    c9 = dblarr(72)
    c10 = dblarr(72)

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

                        readcol, atmname, c1_in, c2_in, c3_in, c4_in, c5_in, c6_in, c7_in, c8_in, c9_in, c10_in, format='D,D,D,D,D,D,D,D,D,D', numline=72, skipline=23, /silent
                        c1 += factor*c1_in
                        c2 += factor*c2_in
                        c3 += factor*c3_in
                        c4 += factor*c4_in
                        c5 += factor*c5_in
                        c6 += factor*c6_in
                        c7 += factor*c7_in
                        c8 += factor*c8_in
                        c9 += factor*c9_in
                        c10 += factor*c10_in
                    endfor
                endfor
            endfor
        endfor
    endfor

    rhox = c1
    mach = c9/c10
    plot, rhox, mach
end
