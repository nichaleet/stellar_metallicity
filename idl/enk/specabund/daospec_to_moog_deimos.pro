pro daospec_to_moog_deimos, mask
    f = file_search(getenv('ahome')+'deimos/'+mask+'_daospec/*.daospec', count=c)
    star = strarr(c)
    specsegment = strarr(c)
    for i=0,c-1 do begin
        fstr = strsplit(file_basename(f[i]), '_.', /extract)
        star[i] = fstr[0]
        specsegment[i] = fstr[1]
    endfor
    stars = star[uniq(star, sort(star))]
    nstars = n_elements(stars)

    for i=0,nstars-1 do begin
        wstar = where(star eq stars[i], c)
        instring = ' '
        moog = {lambda:0d, species:0d, ep:0d, loggf:0d, ew:0d, ewerr:0d, vr:0d}
        moog = replicate(moog, 2000) 
        k = 0L
        for j=0,c-1 do begin
            if (file_info(f[wstar[j]])).size lt 1000 then continue
            openr, lun, f[wstar[j]], /get_lun
            skip_lun, lun, 2, /lines
            while ~eof(lun) do begin
                readf, lun, instring, format='(A114)'
                wave = strmid(instring, 44, 8)
                if strtrim(wave, 2) eq '' then continue
                reads, strmid(instring, 56, 29), species, ep, loggf
                waveobs = double(strmid(instring, 0, 8))
                ew = double(strmid(instring, 19, 8))
                ewerr = double(strmid(instring, 28, 6))
                q = double(strmid(instring, 35, 7))
                if q gt 2 then continue
                moog[k].lambda = wave
                moog[k].species = double(species)
                moog[k].ep = double(ep)
                moog[k].loggf = double(loggf)
                moog[k].ew = ew
                moog[k].ewerr = ewerr
                moog[k].vr = ((waveobs-wave)/ wave) * 2.99792458d5
                k++
            endwhile
            close, lun
            free_lun, lun
        endfor
        if k lt 10 then continue
        moog = moog[0:k-1]

        for j=0,19 do begin
            w = where(abs(moog.vr - mean(moog.vr)) lt 3*stddev(moog.vr), c)
            moog = moog[w]
        endfor
        vr_avg = mean(moog.vr)
        vrerr_avg = stddev(moog.vr) / sqrt(double(c))

        j = -1L
        done = 0
        while 1 do begin
            j++
            if j ge n_elements(moog) then break
            w = where(moog.lambda eq moog[j].lambda and moog.species eq moog[j].species and moog.ep eq moog[j].ep and moog.loggf eq moog[j].loggf, c)
            if c eq 1 then continue
            if c gt 2 then begin
                message, 'I found the same line on more than two orders.', /info
                continue
            endif
            newew = total(moog[w].ew/moog[w].ewerr^2d) / total(1d / moog[w].ewerr^2d)
            newewerr = sqrt(1d / total(1d / moog[w].ewerr^2d))
            moog[j].ew = newew
            moog[j].ewerr = newewerr
            w = w[where(w ne j)]
            ww = complement(w, n_elements(moog))
            moog = moog[ww]
        endwhile
        
        moog = moog[sort(moog.lambda)]
        moog = moog[sort(moog.species)]

        openw, lunw, stars[i]+'.ew', /get_lun
        printf, lunw, stars[i]+string(vr_avg, vrerr_avg, format='("vr=",D+6.1,"+\-",D5.1)')
        for j=0,n_elements(moog)-1 do begin
            if moog[j].ew lt moog[j].ewerr then continue
            if moog[j].ew lt 10 then continue
            printf, lunw, moog[j].lambda, moog[j].species, moog[j].ep, moog[j].loggf, moog[j].ew, format='(D10.3,D10.1,D10.3,D10.3,20(" "),D10.1)'
        endfor
        close, lunw
        free_lun, lunw
    endfor
end
