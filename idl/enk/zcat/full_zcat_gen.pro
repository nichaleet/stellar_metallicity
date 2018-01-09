;+
; NAME:
;       full_zcat_gen
;
; PURPOSE:
;       Generate a zcat with all photometry information, including
;       spatially corrected magnitudes and magnitude errors.	
;
; CALLING SEQUENCE:
;	full_zcat_gen
;
; INPUT:
;	none
;
; OUTPUTS:
;	full_zcat.fits
;
; REVISION HISTORY:
;       enk 2005-11-21
;-

pro full_zcat_gen
    zz = read_zcat()

    pp11 = getphot('11')
    pp12 = getphot('12')
    pp13 = getphot('13')
    pp14 = getphot('14')
    pp21 = getphot('21')
    pp22 = getphot('22')
    pp23 = getphot('23')
    pp31 = getphot('31')
    pp32 = getphot('32')
    pp33 = getphot('33')
    pp41 = getphot('41')
    pp42 = getphot('42')
    pp43 = getphot('43')
    ppinc = [pp14, pp23]
    ppreformed = {objno:ppinc[0].objno, $
                  xs:ppinc[0].xs, $
                  ra:ppinc[0].ra, $
                  dec:ppinc[0].dec, $
                  magu:ppinc[0].magu, $
                  magb:ppinc[0].magb, $
                  magv:ppinc[0].magv, $
                  magr:ppinc[0].magr, $
                  magi:ppinc[0].magi, $
                  magbc:-100.0, $
                  magrc:-100.0, $
                  magic:-100.0, $
                  magb3rg:ppinc[0].magb3rg, $
                  magb6rg:ppinc[0].magb6rg, $
                  magb1s:ppinc[0].magb1s, $
                  magb2s:ppinc[0].magb2s, $
                  magr3rg:ppinc[0].magr3rg, $
                  magr6rg:ppinc[0].magr6rg, $
                  magr1s:ppinc[0].magr1s, $
                  magr2s:ppinc[0].magr2s, $
                  magi3rg:ppinc[0].magi3rg, $
                  magi6rg:ppinc[0].magi6rg, $
                  magi1s:ppinc[0].magi1s, $
                  magi2s:ppinc[0].magi2s, $
                  magberr:ppinc[0].magberr, $
                  magrerr:ppinc[0].magrerr, $
                  magierr:ppinc[0].magierr, $
                  rg:ppinc[0].rg, $
                  nu:ppinc[0].nu, $
                  rh:ppinc[0].rh, $
                  rp:ppinc[0].rp, $
                  rql:ppinc[0].rql, $
                  rqu:ppinc[0].rqu, $
                  nbad:ppinc[0].nbad, $
                  fmax:ppinc[0].fmax, $
                  pa:ppinc[0].pa, $
                  e2:ppinc[0].e2, $
                  quality:ppinc[0].quality, $
                  pgal:ppinc[0].pgal, $
                  photoz:ppinc[0].photoz, $
                  probgt07:ppinc[0].probgt07, $
                  badflag:ppinc[0].badflag, $
                  sfd_ebv:ppinc[0].sfd_ebv}
    ppreformed = replicate(ppreformed, n_elements(ppinc))
    for i=1L, n_elements(ppinc)-1 do begin
        ppreformed[i] = {objno:ppinc[i].objno, $
                      xs:ppinc[i].xs, $
                      ra:ppinc[i].ra, $
                      dec:ppinc[i].dec, $
                      magu:ppinc[i].magu, $
                      magb:ppinc[i].magb, $
                      magv:ppinc[i].magv, $
                      magr:ppinc[i].magr, $
                      magi:ppinc[i].magi, $
                      magbc:-100.0, $
                      magrc:-100.0, $
                      magic:-100.0, $
                      magb3rg:ppinc[i].magb3rg, $
                      magb6rg:ppinc[i].magb6rg, $
                      magb1s:ppinc[i].magb1s, $
                      magb2s:ppinc[i].magb2s, $
                      magr3rg:ppinc[i].magr3rg, $
                      magr6rg:ppinc[i].magr6rg, $
                      magr1s:ppinc[i].magr1s, $
                      magr2s:ppinc[i].magr2s, $
                      magi3rg:ppinc[i].magi3rg, $
                      magi6rg:ppinc[i].magi6rg, $
                      magi1s:ppinc[i].magi1s, $
                      magi2s:ppinc[i].magi2s, $
                      magberr:ppinc[i].magberr, $
                      magrerr:ppinc[i].magrerr, $
                      magierr:ppinc[i].magierr, $
                      rg:ppinc[i].rg, $
                      nu:ppinc[i].nu, $
                      rh:ppinc[i].rh, $
                      rp:ppinc[i].rp, $
                      rql:ppinc[i].rql, $
                      rqu:ppinc[i].rqu, $
                      nbad:ppinc[i].nbad, $
                      fmax:ppinc[i].fmax, $
                      pa:ppinc[i].pa, $
                      e2:ppinc[i].e2, $
                      quality:ppinc[i].quality, $
                      pgal:ppinc[i].pgal, $
                      photoz:ppinc[i].photoz, $
                      probgt07:ppinc[i].probgt07, $
                      badflag:ppinc[i].badflag, $
                      sfd_ebv:ppinc[i].sfd_ebv}
    endfor

    pp = [pp11, pp12, pp13, $
          pp21, pp22, $
          pp31, pp32, pp33, $
          pp41, pp42, pp43, $
          ppreformed]
    match, long(zz.objname), pp.objno, zzw, ppw
    zz = zz[zzw]
    pp = pp[ppw]
    nobj = n_elements(zzw)
    print, 'number of objects in full zcat: ', nobj
    zcat = {zcat, $
            objno:pp[0].objno, $
            ra:pp[0].ra, $
            dec:pp[0].dec, $
            magb:pp[0].magb, $
            magr:pp[0].magr, $
            magi:pp[0].magi, $
            magbc:pp[0].magbc, $
            magrc:pp[0].magrc, $
            magic:pp[0].magic, $
            magberr:pp[0].magberr, $
            magrerr:pp[0].magrerr, $
            magierr:pp[0].magierr, $
            rg:pp[0].rg, $
            rh:pp[0].rh, $
            pgal:pp[0].pgal, $
            photoz:pp[0].photoz, $
            probgt07:pp[0].probgt07, $
            sfd_ebv:pp[0].sfd_ebv, $
            class:zz[0].class, $
            subclass:zz[0].subclass, $
            objname:zz[0].objname, $
            maskname:zz[0].maskname, $
            slitname:zz[0].slitname, $
            date:zz[0].date, $
            mjd:double(zz[0].mjd), $
            z:zz[0].z, $
            z_err:zz[0].z_err, $
            rchi2:zz[0].rchi2, $
            dof:zz[0].dof, $
            rchi2diff:zz[0].rchi2diff, $
            tfile:zz[0].tfile, $
            tcolumn:zz[0].tcolumn, $
            npoly:zz[0].npoly, $
            theta:zz[0].theta, $
            vdisp:zz[0].vdisp, $
            vdisp_err:zz[0].vdisp_err, $
            zquality:zz[0].zquality, $
            comment:zz[0].comment, $
            zspecer:'unknown            '}
    zcat = replicate(zcat, nobj)
    for i=1L, nobj-1 do begin
        zcat[i] = {zcat, $
                   objno:pp[i].objno, $
                   ra:pp[i].ra, $
                   dec:pp[i].dec, $
                   magb:pp[i].magb, $
                   magr:pp[i].magr, $
                   magi:pp[i].magi, $
                   magbc:pp[i].magbc, $
                   magrc:pp[i].magrc, $
                   magic:pp[i].magic, $
                   magberr:pp[i].magberr, $
                   magrerr:pp[i].magrerr, $
                   magierr:pp[i].magierr, $
                   rg:pp[i].rg, $
                   rh:pp[i].rh, $
                   pgal:pp[i].pgal, $
                   photoz:pp[i].photoz, $
                   probgt07:pp[i].probgt07, $
                   sfd_ebv:pp[i].sfd_ebv, $
                   class:zz[i].class, $
                   subclass:zz[i].subclass, $
                   objname:zz[i].objname, $
                   maskname:zz[i].maskname, $
                   slitname:zz[i].slitname, $
                   date:zz[i].date, $
                   mjd:double(zz[i].mjd), $
                   z:zz[i].z, $
                   z_err:zz[i].z_err, $
                   rchi2:zz[i].rchi2, $
                   dof:zz[i].dof, $
                   rchi2diff:zz[i].rchi2diff, $
                   tfile:zz[i].tfile, $
                   tcolumn:zz[i].tcolumn, $
                   npoly:zz[i].npoly, $
                   theta:zz[i].theta, $
                   vdisp:zz[i].vdisp, $
                   vdisp_err:zz[i].vdisp_err, $
                   zquality:zz[i].zquality, $
                   comment:zz[i].comment, $
                   zspecer:'unknown            '}
    endfor
    zcat = zcat[sort(zcat.objno)]
    mwrfits, zcat, 'full_zcat.fits'
end
