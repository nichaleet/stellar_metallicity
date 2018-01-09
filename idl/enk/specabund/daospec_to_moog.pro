pro daospec_to_moog, star, moog=moog, jlc=jlc, seed=seed, mctrial=mctrial
    if keyword_set(jlc) then readcol, getenv('CALTECH')+'linelist/linelist.jlc', lambdaj, speciesj, epj, loggfj, dampingj, skipline=1, format='D,D,D,D,D', /silent
    if keyword_set(mctrial) then readcol, getenv('CALTECH')+'hires/'+star+'/'+star+'.ew', lambdam, speciesm, epm, loggfm, skipline=1, format='D,D,D,D', /silent

    hdr = headfits(getenv('chome')+'hires/'+star+'/'+(keyword_set(mctrial) ? string(mctrial, format='(I05)')+'/' : '')+star+'_030.fits', /silent)
    ras = sxpar(hdr, 'RA')
    decs = sxpar(hdr, 'DEC')
    jd = sxpar(hdr, 'MJD')
    get_coords, coords, instring=ras+' '+decs
    ra = coords[0]*15.
    dec = coords[1]
    heliocorr = helio_deimos(ra, dec, 2000.0, jd=double(jd) + 2400000.5d)
    applyhelio = star eq 'S150' or star eq 'S24' or star eq 'S53' or star eq 'S7'

    f = file_search(getenv('chome')+'hires/'+star+'/'+(keyword_set(mctrial) ? string(mctrial, format='(I05)')+'/' : '')+star+'_*.daospec', count=c)
    instring = ' '
    daodata = replicate({lambda:0d, ew:0d, ewerr:0d, order:0}, 100000)
    moog = {lambda:0d, species:0d, ep:0d, loggf:0d, damping:2d, ew:0d, ewerr:0d, vr:0d, order:0}
    moog = replicate(moog, 10000)
    j = 0L
    k = 0L
    for i=0,c-1 do begin
        if (file_info(f[i])).size lt 100 then continue
        openr, lun, f[i], /get_lun
        skip_lun, lun, 2, /lines
        while ~eof(lun) do begin
            readf, lun, instring, format='(A114)'
            daowave = strmid(instring, 10, 8)
            wave = strmid(instring, 44, 8)
            ew = double(strmid(instring, 19, 8))
            ewerr = double(strmid(instring, 28, 6))
            quality = double(strmid(instring, 35, 7))
            if quality gt 20 then continue
            daodata[k].lambda = daowave
            daodata[k].ew = ew
            daodata[k].ewerr = ewerr
            daodata[k].order = i+1
            k++
            if strtrim(wave, 2) eq '' then continue
            reads, strmid(instring, keyword_set(moog) ? 56 : 53, 29), species, ep, loggf
            waveobs = double(strmid(instring, 0, 8))
            moog[j].lambda = wave
            moog[j].species = double(species)
            moog[j].ep = double(ep)
            moog[j].loggf = double(loggf)
            moog[j].ew = ew
            moog[j].ewerr = ewerr
            moog[j].vr = ((waveobs-wave)/ wave) * 2.99792458d5 - (applyhelio ? heliocorr : 0d)
            moog[j].order = i+1
            j++
        endwhile
        close, lun
        free_lun, lun
    endfor
    moog = moog[0:j-1]
    daodata = daodata[0:k-1]
    moog = moog[where(moog.ewerr gt 0)]
    daodata = daodata[where(daodata.ewerr gt 0)]

    for i=0,19 do begin
        if i lt 5 then w = where(abs(moog.vr - mean(moog.vr)) lt 3*stddev(moog.vr), c) else w = where(abs(moog.vr - mean(moog.vr)) lt 3*stddev(moog.vr) and abs(moog.vr - median(moog.vr)) lt 10, c)
        moog = moog[w]
    endfor
    vr_avg = mean(moog.vr)
    vrerr_avg = stddev(moog.vr) / sqrt(double(c))
    if vrerr_avg gt 1 then message, 'Bad line identifications found.'
    vr = vr_avg
    vrerr = vrerr_avg
    save, vr, vrerr, filename=getenv('CALTECH')+'hires/'+star+'_vr.sav'

    ;openw, lun2, 'n5024_vr.dat', /get_lun, /append
    ;printf, lun2, star, vr_avg, vrerr_avg, ras, decs, format='(A-15,1X,D+8.3,1X,D7.3,1X,A15,1X,A15)'
    ;close, lun2
    ;free_lun, lun2

    if keyword_set(jlc) then begin
        keep = bytarr(j)
        for i=0,c-1 do begin
            w = where(abs(lambdaj-moog[i].lambda) lt 0.03 and abs(speciesj-moog[i].species) lt 0.01 and abs(epj-moog[i].ep) lt 0.1, c)
            if c eq 1 then begin
                keep[i] = 1
                moog[i].damping = dampingj[w]
            endif
        endfor

        w = where(keep, c)
        if c eq 0 then message, 'No lines in common.'
        moog = moog[w]
    endif

    if keyword_set(mctrial) then begin
        keep = bytarr(c)
        for i=0,c-1 do begin
            w = where(abs(lambdam-moog[i].lambda) lt 0.02 and abs(speciesm-moog[i].species) lt 0.01 and abs(epm-moog[i].ep) lt 0.1, c)
            if c eq 1 then begin
                keep[i] = 1
            endif
        endfor

        w = where(keep, c)
        if c eq 0 then message, 'No lines in common.'
        moog = moog[w]
    endif

    for i=0,c-1 do begin
        w = where(abs(daodata.lambda-moog[i].lambda) lt 0.05 and daodata.order eq moog[i].order, cclose)
        if cclose gt 1 then begin
            ;print, moog[i].lambda, moog[i].species, moog[i].ep, moog[i].loggf, daodata[w].ew, format='(D10.3,D10.1,D10.3,D10.3,10(D6.1))'
            moog[i].ew = total(daodata[w].ew)
            moog[i].ewerr = sqrt(total((daodata[w].ewerr)^2d))
        endif
    endfor

    i = -1L
    done = 0
    while 1 do begin
        i++
        if i ge n_elements(moog) then break
        w = where(moog.lambda eq moog[i].lambda and moog.species eq moog[i].species and moog.ep eq moog[i].ep and moog.loggf eq moog[i].loggf, c)
        if c eq 1 then continue
        if c gt 2 then message, 'I found the same line on more than two orders.'
        newew = total(moog[w].ew/moog[w].ewerr^2d) / total(1d / moog[w].ewerr^2d)
        newewerr = sqrt(1d / total(1d / moog[w].ewerr^2d))
        moog[i].ew = newew
        moog[i].ewerr = newewerr
        w = w[where(w ne i)]
        ww = complement(w, n_elements(moog))
        moog = moog[ww]
    endwhile
    
    moog = moog[where(moog.ew gt 10 and moog.ew gt 3.*moog.ewerr)]
    moog = moog[sort(moog.lambda)]
    moog = moog[sort(moog.species)]

    case 1 of
        keyword_set(seed): suffix =  '_mc'
        keyword_set(mctrial): suffix = '_'+string(mctrial, format='(I05)')
        else: suffix = ''
    endcase
    openw, lunw, getenv('CALTECH')+'hires/'+star+'/'+star+suffix+'.ew', /get_lun
    printf, lunw, star
    for i=0,n_elements(moog)-1 do begin
        if moog[i].ew lt 3*moog[i].ewerr then continue
        if moog[i].ew lt 10 then continue
        ew_i = keyword_set(seed) ? moog[i].ew+moog[i].ewerr*randomn(seed, 1) : moog[i].ew
        printf, lunw, moog[i].lambda, moog[i].species, moog[i].ep, moog[i].loggf, moog[i].damping, ew_i, format='(D10.3,D10.1,D10.3,D10.3,G10.4,10(" "),D10.1)'
    endfor
    close, lunw
    free_lun, lunw
end
