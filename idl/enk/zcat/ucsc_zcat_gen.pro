;+
; NAME:
;       
;
; PURPOSE:
;	It is difficult to find an updated zcat for those masks
;	zspec'ed at Santa Cruz.  This routine generates a structure
;	called ucsc_zcat, which is a Santa Cruz zcat.  It will always be
;	current as long as no directory is added to the list of
;	directories which contain Santa Cruz zspec results.
;       
; CALLING SEQUENCE:
;	ucsc_zcat_gen
;
; INPUT:
;	none
;
; OUTPUTS:
;	ucsc_zcat.fits
;
; REVISION HISTORY:
;       enk 2005-11-11
;-

pro ucsc_zcat_gen
    path = strarr(39)
    name = strarr(39)
    path[0] = '/net/asturias/b/kassin/zspec_results/anne/'
    name[0] = 'Anne Metevier'
    path[1] = '/net/asturias/b/kassin/zspec_results/bjw/'
    name[1] = 'Benjamin Weiner'
    path[2] = '/net/asturias/b/kassin/zspec_results/cnaw/'
    name[2] = 'Christopher Willmer'
    path[3] = '/net/asturias/b/kassin/zspec_results/cnaw.gst131/'
    name[3] = 'Christopher Willmer'
    path[4] = '/net/asturias/b/kassin/zspec_results/kai/'
    name[4] = 'Kai Noeske'
    path[5] = '/net/asturias/b/kassin/zspec_results/koo/'
    name[5] = 'David Koo'
    path[6] = '/net/asturias/b/kassin/zspec_results/lemaux/'
    name[6] = 'Brian Lemaux'
    path[7] = '/net/asturias/b/kassin/zspec_results/lihwai/'
    name[7] = 'Lihwai Lin'
    path[8] = '/net/asturias/b/kassin/zspec_results/npk/'
    name[8] = 'Nick Konidaris'
    path[9] = '/net/asturias/a/kassin/zspec04/gen/'
    name[9] = 'Genevieve Graves'
    path[10] = '/net/asturias/a/kassin/zspec04/justin/'
    name[10] = 'Justin Harker'
    path[11] = '/net/asturias/a/kassin/zspec04/kai/'
    name[11] = 'Kai Noeske'
    path[12] = '/net/asturias/a/kassin/zspec04/lihwai/'
    name[12] = 'Lihwai Lin'
    path[13] = '/net/asturias/a/kassin/zspec04/npk/'
    name[13] = 'Nick Konidaris'
    path[14] = '/net/asturias/a/kassin/zspec05/anne/'
    name[14] = 'Anne Metevier'
    path[15] = '/net/asturias/a/kassin/zspec05/cnaw/'
    name[15] = 'Christopher Willmer'
    path[16] = '/net/asturias/a/kassin/zspec05/evan/'
    name[16] = 'Evan Kirby'
    path[17] = '/net/asturias/a/kassin/zspec05/faber/'
    name[17] = 'Sandy Faber'
    path[18] = '/net/asturias/a/kassin/zspec05/faber_prev05/'
    name[18] = 'Sandy Faber'
    path[19] = '/net/asturias/a/kassin/zspec05/gen/'
    name[19] = 'Genevieve Graves'
    path[20] = '/net/asturias/a/kassin/zspec05/jmel/'
    name[20] = 'Jason Melbourne'
    path[21] = '/net/asturias/a/kassin/zspec05/justin/'
    name[21] = 'Justin Harker'
    path[22] = '/net/asturias/a/kassin/zspec05/kai/'
    name[22] = 'Kai Noeske'
    path[23] = '/net/asturias/a/kassin/zspec05/kassin/'
    name[23] = 'Susan Kassin'
    path[24] = '/net/asturias/a/kassin/zspec05/koo/'
    name[24] = 'David Koo'
    path[25] = '/net/asturias/a/kassin/zspec05/lihwai/'
    name[25] = 'Lihwai Lin'
    path[26] = '/net/asturias/a/kassin/zspec05/lotz/'
    name[26] = 'Jennifer Lotz'
    path[27] = '/net/asturias/a/kassin/zspec05/npk/'
    name[27] = 'Nick Konidaris'
    path[28] = '/net/asturias/a/kassin/zspec05/raja/'
    name[28] = 'Raja Guha Thakurta'
    path[29] = '/net/asturias/a/bjw/zspecaug03/bjw/'
    name[29] = 'Benjamin Weiner'
    path[30] = '/net/asturias/a/bjw/zspecaug03/brian/'
    name[30] = 'Brian Lemaux'
    path[31] = '/net/asturias/a/bjw/zspecaug03/cnaw/'
    name[31] = 'Christopher Willmer'
    path[32] = '/net/asturias/a/bjw/zspecaug03/deepteam/'
    name[32] = 'Brian Lemaux'
    path[33] = '/net/asturias/a/bjw/zspecaug03/gst131/'
    name[33] = 'gst131'
    path[34] = '/net/asturias/a/bjw/zspecaug03/lihwai/'
    name[34] = 'Lihwai Lin'
    path[35] = '/net/asturias/a/bjw/zspecaug03/npk/'
    name[35] = 'Nick Konidaris'
    path[36] = '/net/asturias/a/bjw/zspecaug03/raja/'
    name[36] = 'Raja Guha Thakurta'
    path[37] = '/net/pemla/u/raja/DEEP/zspec/egs_acs/outp_fits_sav/set1/'
    name[37] = 'Raja Guha Thakurta'
    path[38] = '/net/pemla/u/raja/DEEP/zspec/egs_acs/outp_fits_sav/set2/'
    name[38] = 'Raja Guha Thakurta'

    zzcreated = 0
    for i=0L,n_elements(path)-1 do begin
        files = findfile(path[i]+'*.sav')
        print, path[i]
        for j=0L,n_elements(files)-1 do begin
            ;print, 'i: ', i, '  j: ', j
            restore, files[j]
            if zzcreated then begin
                zz = [zz, results]
                zspecer = [zspecer, replicate(name[i], n_elements(results))]
            endif else begin
                zz = results
                zspecer = replicate(name[i], n_elements(results))
                zzcreated = 1
            endelse
        endfor
    endfor

    ;mn = strcompress(zz.maskname, /remove_all)
    ;maskabbr = strmid(mn, 0, 2)
    ;maskabbr = maskabbr[sort(maskabbr)]
    ;maskabbruniq = maskabbr[uniq(maskabbr)]
    ;npp = n_elements(maskabbruniq)
    ;for i=0, npp-1 do begin
    ;    print, maskabbruniq[i]
    ;    if i eq 0 then begin
    ;        pp = getphot(maskabbruniq[i])
    ;    endif else begin
    ;        pp = [pp, getphot(maskabbruniq[i])]
    ;    endelse            
    ;endfor
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
    for i=long64(1), n_elements(ppinc)-1 do begin
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
    zspecer = zspecer[zzw]
    pp = pp[ppw]
    nobj = n_elements(zzw)
    print, 'number of objects in UCSC zcat: ', nobj
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
            class:zz[0].zresult.class, $
            subclass:zz[0].zresult.subclass, $
            objname:zz[0].objname, $
            maskname:zz[0].maskname, $
            slitname:zz[0].slitname, $
            date:zz[0].zresult.date, $
            mjd:zz[0].zresult.mjd, $
            z:zz[0].zresult.z, $
            z_err:zz[0].zresult.z_err, $
            rchi2:zz[0].zresult.rchi2, $
            dof:zz[0].zresult.dof, $
            rchi2diff:zz[0].zresult.rchi2diff, $
            tfile:zz[0].zresult.tfile, $
            tcolumn:zz[0].zresult.tcolumn, $
            npoly:zz[0].zresult.npoly, $
            theta:zz[0].zresult.theta, $
            vdisp:zz[0].zresult.vdisp, $
            vdisp_err:zz[0].zresult.vdisp_err, $
            zquality:zz[0].zresult.zquality, $
            comment:zz[0].zresult.comment, $
            zspecer:zspecer[0]}            
    zcat = replicate(zcat, nobj)
    for i=long64(1), nobj-1 do begin
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
                   class:zz[i].zresult.class, $
                   subclass:zz[i].zresult.subclass, $
                   objname:zz[i].objname, $
                   maskname:zz[i].maskname, $
                   slitname:zz[i].slitname, $
                   date:zz[i].zresult.date, $
                   mjd:zz[i].zresult.mjd, $
                   z:zz[i].zresult.z, $
                   z_err:zz[i].zresult.z_err, $
                   rchi2:zz[i].zresult.rchi2, $
                   dof:zz[i].zresult.dof, $
                   rchi2diff:zz[i].zresult.rchi2diff, $
                   tfile:zz[i].zresult.tfile, $
                   tcolumn:zz[i].zresult.tcolumn, $
                   npoly:zz[i].zresult.npoly, $
                   theta:zz[i].zresult.theta, $
                   vdisp:zz[i].zresult.vdisp, $
                   vdisp_err:zz[i].zresult.vdisp_err, $
                   zquality:zz[i].zresult.zquality, $
                   comment:zz[i].zresult.comment, $
                   zspecer:zspecer[i]}
    endfor
    zcat = zcat[sort(zcat.objno)]
    ;save, zcat, filename='ucsc_zcat.sav'
    mwrfits, zcat, 'ucsc_zcat.fits'
end
