pro kcorr1
    zcat = all_zcat()
    openw, 1, 'zcat.in'
    for i=0L, n_elements(zcat)-1 do printf, 1, zcat[i].objno, zcat[i].z, zcat[i].zquality, zcat[i].magB, zcat[i].magR, zcat[i].magI
    close, 1
end

pro kcorr2
    readcol, 'zcat.out', objno, z, magR, BR, RI, K_RB, MB0, UB0, UV0, BV0, BR0, UG0, GR0, u2000b, u2500b, u2800b, UB_interp, B_interp, BV_interp, mass, format='L,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D'
    zcat = all_zcat()
    zcat = zcat[uniq(zcat.objno, sort(zcat.objno))]
    w = uniq(objno, sort(objno))
    objno = objno[w]
    K_RB = K_RB[w]
    MB0 = MB0[w]
    UB0 = UB0[w]
    UV0 = UV0[w]
    BV0 = BV0[w]
    BR0 = BR0[w]
    UG0 = UG0[w]
    GR0 = GR0[w]
    mass = mass[w]
    match, zcat.objno, objno, w1, w2
    j = 0L
    if contains(w1, 0) then begin
        kzcat = {kzcat, $
                 objno:zcat[0].objno, $
                 ra:zcat[0].ra, $
                 dec:zcat[0].dec, $
                 magb:zcat[0].magb, $
                 magr:zcat[0].magr, $
                 magi:zcat[0].magi, $
                 magbc:zcat[0].magbc, $
                 magrc:zcat[0].magrc, $
                 magic:zcat[0].magic, $
                 magberr:zcat[0].magberr, $
                 magrerr:zcat[0].magrerr, $
                 magierr:zcat[0].magierr, $
                 rg:zcat[0].rg, $
                 rh:zcat[0].rh, $
                 pgal:zcat[0].pgal, $
                 photoz:zcat[0].photoz, $
                 probgt07:zcat[0].probgt07, $
                 sfd_ebv:zcat[0].sfd_ebv, $
                 class:zcat[0].class, $
                 subclass:zcat[0].subclass, $
                 objname:zcat[0].objname, $
                 maskname:zcat[0].maskname, $
                 slitname:zcat[0].slitname, $
                 date:zcat[0].date, $
                 mjd:zcat[0].mjd, $
                 z:zcat[0].z, $
                 z_err:zcat[0].z_err, $
                 rchi2:zcat[0].rchi2, $
                 dof:zcat[0].dof, $
                 rchi2diff:zcat[0].rchi2diff, $
                 tfile:zcat[0].tfile, $
                 tcolumn:zcat[0].tcolumn, $
                 npoly:zcat[0].npoly, $
                 theta:zcat[0].theta, $
                 vdisp:zcat[0].vdisp, $
                 vdisp_err:zcat[0].vdisp_err, $
                 zquality:zcat[0].zquality, $
                 comment:zcat[0].comment, $
                 zspecer:zcat[0].zspecer, $
                 K_RB:K_RB[w2[j]], $
                 MB:MB0[w2[j]], $
                 UB:UB0[w2[j]], $
                 UV:UV0[w2[j]], $
                 BV:BV0[w2[j]], $
                 BR:BR0[w2[j]], $
                 UG:UG0[w2[j]], $
                 GR:GR0[w2[j]], $
                 stelmass:mass[w2[j]]}
        j = j + 1
    endif else begin
        kzcat = {kzcat, $
                 objno:zcat[0].objno, $
                 ra:zcat[0].ra, $
                 dec:zcat[0].dec, $
                 magb:zcat[0].magb, $
                 magr:zcat[0].magr, $
                 magi:zcat[0].magi, $
                 magbc:zcat[0].magbc, $
                 magrc:zcat[0].magrc, $
                 magic:zcat[0].magic, $
                 magberr:zcat[0].magberr, $
                 magrerr:zcat[0].magrerr, $
                 magierr:zcat[0].magierr, $
                 rg:zcat[0].rg, $
                 rh:zcat[0].rh, $
                 pgal:zcat[0].pgal, $
                 photoz:zcat[0].photoz, $
                 probgt07:zcat[0].probgt07, $
                 sfd_ebv:zcat[0].sfd_ebv, $
                 class:zcat[0].class, $
                 subclass:zcat[0].subclass, $
                 objname:zcat[0].objname, $
                 maskname:zcat[0].maskname, $
                 slitname:zcat[0].slitname, $
                 date:zcat[0].date, $
                 mjd:zcat[0].mjd, $
                 z:zcat[0].z, $
                 z_err:zcat[0].z_err, $
                 rchi2:zcat[0].rchi2, $
                 dof:zcat[0].dof, $
                 rchi2diff:zcat[0].rchi2diff, $
                 tfile:zcat[0].tfile, $
                 tcolumn:zcat[0].tcolumn, $
                 npoly:zcat[0].npoly, $
                 theta:zcat[0].theta, $
                 vdisp:zcat[0].vdisp, $
                 vdisp_err:zcat[0].vdisp_err, $
                 zquality:zcat[0].zquality, $
                 comment:zcat[0].comment, $
                 zspecer:zcat[0].zspecer, $
                 K_RB:-100d, $
                 MB:-100d, $
                 UB:-100d, $
                 UV:-100d, $
                 BV:-100d, $
                 BR:-100d, $
                 UG:-100d, $
                 GR:-100d, $
                 stelmass:-100d}
    endelse
    kzcat = replicate(kzcat, n_elements(zcat))
    for i=1L, n_elements(zcat)-1 do begin
        if contains(w1, i) then begin
            kzcat[i] = {kzcat, $
                     objno:zcat[i].objno, $
                     ra:zcat[i].ra, $
                     dec:zcat[i].dec, $
                     magb:zcat[i].magb, $
                     magr:zcat[i].magr, $
                     magi:zcat[i].magi, $
                     magbc:zcat[i].magbc, $
                     magrc:zcat[i].magrc, $
                     magic:zcat[i].magic, $
                     magberr:zcat[i].magberr, $
                     magrerr:zcat[i].magrerr, $
                     magierr:zcat[i].magierr, $
                     rg:zcat[i].rg, $
                     rh:zcat[i].rh, $
                     pgal:zcat[i].pgal, $
                     photoz:zcat[i].photoz, $
                     probgt07:zcat[i].probgt07, $
                     sfd_ebv:zcat[i].sfd_ebv, $
                     class:zcat[i].class, $
                     subclass:zcat[i].subclass, $
                     objname:zcat[i].objname, $
                     maskname:zcat[i].maskname, $
                     slitname:zcat[i].slitname, $
                     date:zcat[i].date, $
                     mjd:zcat[i].mjd, $
                     z:zcat[i].z, $
                     z_err:zcat[i].z_err, $
                     rchi2:zcat[i].rchi2, $
                     dof:zcat[i].dof, $
                     rchi2diff:zcat[i].rchi2diff, $
                     tfile:zcat[i].tfile, $
                     tcolumn:zcat[i].tcolumn, $
                     npoly:zcat[i].npoly, $
                     theta:zcat[i].theta, $
                     vdisp:zcat[i].vdisp, $
                     vdisp_err:zcat[i].vdisp_err, $
                     zquality:zcat[i].zquality, $
                     comment:zcat[i].comment, $
                     zspecer:zcat[i].zspecer, $
                     K_RB:K_RB[w2[j]], $
                     MB:MB0[w2[j]], $
                     UB:UB0[w2[j]], $
                     UV:UV0[w2[j]], $
                     BV:BV0[w2[j]], $
                     BR:BR0[w2[j]], $
                     UG:UG0[w2[j]], $
                     GR:GR0[w2[j]], $
                     stelmass:mass[w2[j]]}
            j = j + 1
        endif else begin
            kzcat[i] = {kzcat, $
                     objno:zcat[i].objno, $
                     ra:zcat[i].ra, $
                     dec:zcat[i].dec, $
                     magb:zcat[i].magb, $
                     magr:zcat[i].magr, $
                     magi:zcat[i].magi, $
                     magbc:zcat[i].magbc, $
                     magrc:zcat[i].magrc, $
                     magic:zcat[i].magic, $
                     magberr:zcat[i].magberr, $
                     magrerr:zcat[i].magrerr, $
                     magierr:zcat[i].magierr, $
                     rg:zcat[i].rg, $
                     rh:zcat[i].rh, $
                     pgal:zcat[i].pgal, $
                     photoz:zcat[i].photoz, $
                     probgt07:zcat[i].probgt07, $
                     sfd_ebv:zcat[i].sfd_ebv, $
                     class:zcat[i].class, $
                     subclass:zcat[i].subclass, $
                     objname:zcat[i].objname, $
                     maskname:zcat[i].maskname, $
                     slitname:zcat[i].slitname, $
                     date:zcat[i].date, $
                     mjd:zcat[i].mjd, $
                     z:zcat[i].z, $
                     z_err:zcat[i].z_err, $
                     rchi2:zcat[i].rchi2, $
                     dof:zcat[i].dof, $
                     rchi2diff:zcat[i].rchi2diff, $
                     tfile:zcat[i].tfile, $
                     tcolumn:zcat[i].tcolumn, $
                     npoly:zcat[i].npoly, $
                     theta:zcat[i].theta, $
                     vdisp:zcat[i].vdisp, $
                     vdisp_err:zcat[i].vdisp_err, $
                     zquality:zcat[i].zquality, $
                     comment:zcat[i].comment, $
                     zspecer:zcat[i].zspecer, $
                     K_RB:-100d, $
                     MB:-100d, $
                     UB:-100d, $
                     UV:-100d, $
                     BV:-100d, $
                     BR:-100d, $
                     UG:-100d, $
                     GR:-100d, $
                     stelmass:-100d}
        endelse        
    endfor
    mwrfits, kzcat, 'all_kzcat.fits'
end
