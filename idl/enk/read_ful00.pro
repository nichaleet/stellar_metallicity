function read_ful00
    ful00 = {HIP:0L, HD:0L, name:' ', othername:' ', ra:0d, dec:0d, pmra:0d, pmdec:0d, sptype:' ', v:0d, texp:0d, snr:0d, teff:0L, logg:0d, vt:0d, feh:0d, ali:0d, nafe:0d, mgfe:0d, alfe:0d, sife:0d, cafe:0d, tife:0d, vfe:0d, crfe:0d, nife:0d, yfe:0d, zrfe:0d, bafe:0d, eufe:0d}
    ful00 = replicate(ful00, 1000)

    n = 0
    hip = 0L
    hd = 0L
    name = ' '
    othername = ' '
    openr, lun, getenv('M31')+'papers/ful00/ful00_tab2.dat', /get_lun
    skip_lun, lun, 20, /lines
    while ~eof(lun) do begin
        readf, lun, hip, hd, name, othername, v, texp, snr, format='(I6,7X,I6,1X,A11,1X,A9,1X,D4,1X,D4,1X,D3)'
        w = where(ful00.hip eq hip, c)
        if c eq 0 then begin
            ful00[n].hip = hip
            ful00[n].hd = hd
            ful00[n].name = strtrim(name, 2)
            ful00[n].othername = strtrim(othername, 2)
            ful00[n].v = v
            ful00[n].texp = texp
            ful00[n].snr = snr
            n++
        endif else begin
            ful00[w].texp += texp
            ful00[w].snr = sqrt(ful00[w].snr^2. + snr^2.)
        endelse
    endwhile
    close, lun
    free_lun, lun
    ful00 = ful00[0:n-1]

    teff = 0L
    openr, lun, getenv('M31')+'papers/ful00/ful00_tab5.dat', /get_lun
    skip_lun, lun, 23, /lines
    while ~eof(lun) do begin
        readf, lun, hip, teff, logg, feh, vt, format='(I6,24X,I4,1X,D3,1X,D4,1X,D4)'
        w = where(ful00.hip eq hip, c)
        if c eq 1 then begin
            ful00[w].teff = teff
            ful00[w].logg = logg
            ;ful00[w].feh = feh
            ful00[w].vt = vt
        endif
    endwhile
    close, lun
    free_lun, lun

    hip = ' '
    ali = ' '
    nafe = ' '
    mgfe = ' '
    alfe = ' '
    sife = ' '
    cafe = ' '
    tife = ' '
    vfe = ' '
    crfe = ' '
    nife = ' '
    yfe = ' '
    zrfe = ' '
    bafe = ' '
    eufe = ' '
    openr, lun, getenv('M31')+'papers/ful00/ful00_tab6.dat', /get_lun
    skip_lun, lun, 32, /lines
    while ~eof(lun) do begin
        readf, lun, hip, feh, ali, nafe, mgfe, alfe, sife, cafe, tife, vfe, crfe, nife, yfe, zrfe, bafe, eufe, format='(A7,8(1X,A5),1X,7(1X,A5))'
        w = where(strtrim(ful00.hip, 2) eq strtrim(hip, 2), c)
        if c eq 1 then begin
            ful00[w].feh = feh
            ful00[w].ali = strtrim(ali, 2) eq '' ? -999.0 : ali
            ful00[w].nafe = strtrim(nafe, 2) eq '' ? -999.0 : nafe
            ful00[w].mgfe = strtrim(mgfe, 2) eq '' ? -999.0 : mgfe
            ful00[w].alfe = strtrim(alfe, 2) eq '' ? -999.0 : alfe
            ful00[w].sife = strtrim(sife, 2) eq '' ? -999.0 : sife
            ful00[w].cafe = strtrim(cafe, 2) eq '' ? -999.0 : cafe
            ful00[w].tife = strtrim(tife, 2) eq '' ? -999.0 : tife
            ful00[w].vfe = strtrim(vfe, 2) eq '' ? -999.0 : vfe
            ful00[w].crfe = strtrim(crfe, 2) eq '' ? -999.0 : crfe
            ful00[w].nife = strtrim(nife, 2) eq '' ? -999.0 : nife
            ful00[w].yfe = strtrim(yfe, 2) eq '' ? -999.0 : yfe
            ful00[w].zrfe = strtrim(zrfe, 2) eq '' ? -999.0 : zrfe
            ful00[w].bafe = strtrim(bafe, 2) eq '' ? -999.0 : bafe
            ful00[w].eufe = strtrim(eufe, 2) eq '' ? -999.0 : eufe
        endif
    endwhile
    close, lun
    free_lun, lun

    ra = ' '
    dec = ' '
    sptype = ' '
    openr, lun, getenv('M31')+'papers/ful00/ful00_radec.dat', /get_lun
    while ~eof(lun) do begin
        readf, lun, hip, ra, dec, pmra, pmdec, v, sptype, format='(I6,3X,A14,2X,A14,2X,D8,1X,D8,1X,D6,A15)'
        w = where(strtrim(ful00.hip, 2) eq strtrim(hip, 2), c)
        if c eq 1 then begin
            get_coords, coords, instring=ra+' '+dec
            ful00[w].ra = coords[0]*15d
            ful00[w].dec = coords[1]
            ful00[w].pmra = pmra
            ful00[w].pmdec = pmdec
            ful00[w].v = v
            ful00[w].sptype = strtrim(sptype, 2)
        endif
    endwhile
    close, lun
    free_lun, lun

    return, ful00
end
