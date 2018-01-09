function discrepancy, bq, scq
    ucbzc = read_zcat()
    ucsczc = ucsc_zcat()
    listmatch, ucbzc.objno, ucsczc.objno, w1, w2
    ucbzc = ucbzc[w1]
    ucsczc = ucsczc[w2]
    w = where(ucbzc.zquality eq bq)
    ucbzcqi = ucbzc[w]
    listmatch, ucbzcqi.objno, ucsczc.objno, w1, w2
    ucsczcqi = ucsczc[w2]
    w = where(ucsczcqi.zquality eq scq, cwsc)
    if cwsc gt 0 then return, ucsczcqi[w] else return, -1
end

pro qcompare
    ucbzc = read_zcat()
    ucsczc = ucsc_zcat()
    listmatch, ucbzc.objno, ucsczc.objno, w1, w2
    ucbzc = ucbzc[w1]
    ucsczc = ucsczc[w2]

    scq = intarr(7, 7)
    btot = intarr(7)
    for i=-2,4 do begin
        w = where(ucbzc.zquality eq i, cwb)
        btot[i+2] = cwb
        print, 'Berkeley Q = ' + strcompress(i, /remove_all) + ' (' + strcompress(cwb, /remove_all) + ')'
        if cwb gt 0 then begin
            ucbzcqi = ucbzc[w]
            listmatch, ucbzcqi.objno, ucsczc.objno, w1, w2
            ucsczcqi = ucsczc[w2]
            for j=-2,4 do begin
                w = where(ucsczcqi.zquality eq j, cwsc)
                scq[i+2, j+2] = cwsc
                print, 'Santa Cruz Q = ' + strcompress(j, /remove_all) + ': ' + strcompress(cwsc, /remove_all)
            endfor
        endif
        print, ' '
    endfor

    bq = intarr(7, 7)
    sctot = intarr(7)
    for i=-2,4 do begin
        w = where(ucsczc.zquality eq i, cwsc)
        sctot[i+2] = cwsc
        print, 'Santa Cruz Q = ' + strcompress(i, /remove_all) + ' (' + strcompress(cwsc, /remove_all) + ')'
        if cwsc gt 0 then begin
            ucsczcqi = ucbzc[w]
            listmatch, ucsczcqi.objno, ucbzc.objno, w1, w2
            ucbzcqi = ucbzc[w2]
            for j=-2,4 do begin
                w = where(ucbzcqi.zquality eq j, cwb)
                bq[i+2, j+2] = cwb
                print, 'Berkeley Q = ' + strcompress(j, /remove_all) + ': ' + strcompress(cwb, /remove_all)
            endfor
        endif
        print, ' '
    endfor

    openw, 1, 'qcompare.txt'
    printf, 1, 'Berkeley (across) -> Santa Cruz (down)'
    printf, 1, 'Q       ', -2, -1, 0, 1, 2, 3, 4
    for i=-2, 4 do begin
        printf, 1, i, scq[0, i+2], scq[1, i+2], scq[2, i+2], scq[3, i+2], scq[4, i+2], scq[5, i+2], scq[6, i+2]
    endfor
    printf, 1, 'total   ', btot[0], btot[1], btot[2], btot[3], btot[4], btot[5], btot[6]
    printf, 1, ' '
    printf, 1, 'Santa Cruz (across) -> Berkeley (down)'
    printf, 1, 'Q       ', -2, -1, 0, 1, 2, 3, 4
    for i=-2, 4 do begin
        printf, 1, i, bq[0, i+2], bq[1, i+2], bq[2, i+2], bq[3, i+2], bq[4, i+2], bq[5, i+2], bq[6, i+2]
    endfor
    printf, 1, 'total   ', sctot[0], sctot[1], sctot[2], sctot[3], sctot[4], sctot[5], sctot[6]
    close, 1
end
