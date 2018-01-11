pro consolidate, mask, zspecpath=zspecpath, spec1dpath=spec1dpath, alldatapath=alldatapath, photpath=photpath, savflag=savflag
    print, mask

    ;default zspecpath
    if ~keyword_set(zspecpath) then zspecpath = '/net/pemla/data1/data2/m31_deimos/work/zresult/'
    if strmid(zspecpath, 0, 1, /reverse_offset) ne '/' and ~strmatch(zspecpath, '*.fits') then zspecpath += '/'

    ;default spec1dpath
    if ~keyword_set(spec1dpath) then spec1dpath = '/net/pemla/data1/data2/m31_deimos/work'
    ;if strmid(spec1dpath, 0, 1, /reverse_offset) eq '/' then spec1dpath = strmid(spec1dpath, 0, strlen(spec1dpath)-1)

    ;default alldatapath
    ;if ~keyword_set(alldatapath) then alldatapath = '/net/host/something/'
    ;if ~keyword_set(alldatapath) then alldatapath = '/net/schechter/a/ekirby/DEIMOS/phot/'
    if ~keyword_set(alldatapath) then alldatapath = '/net/pemla/data1/data3/kgilbert/DEIMOS/phot/'
    if strmid(alldatapath, 0, 1, /reverse_offset) ne '/' then alldatapath += '/'

    ;extract the field name from the mask name
    field = (strsplit(mask, '_', /extract))[0]

    ;find photometry file
    if ~keyword_set(photpath) then photpath = '/net/pemla/u/raja/m31_deimos/data_status/redux/f07/work/phot/'
    case 1 of
        strmatch(photpath, '*.phot'): begin
            photfile = photpath
            photformat = 'phot'
        end
        strmatch(photpath, '*.cat'): begin
            photfile = photpath
            photformat = 'cat'
        end
        strmatch(photpath, '*.gprob'): begin
            photfile = photpath
            photformat = 'gprob'
        end
        else: begin
            if strmid(photpath, 0, 1, /reverse_offset) ne '/' then photpath += '/'
            photfile = photpath+'all.'+field+'.prob'
            if ~file_test(photfile) then photfile = photpath+'ALL.'+field+'.prob'
            if ~file_test(photfile) then photfile = photpath+'all.'+field+'.prob.txt'
            if ~file_test(photfile) then message, 'Photometry file not found.'
            photformat = 'prob'
        end
    endcase

    ;read the zspec FITS file
    if strmatch(zspecpath, '*.fits') then zspecfile = zspecpath else zspecfile = zspecpath+'zspec.'+mask+'.fits'
    if ~file_test(zspecfile) then message, 'zspecfile not found'
    zspec = mrdfits(zspecfile, 1, /silent)
    if strmid(mask, 0, 1) eq 'd' then zspec = zspec[where(strtrim(zspec.maskname, 2) eq strtrim(mask, 2))]
    tags = tag_names(zspec)
    wstarid = where(strtrim(tags, 2) eq 'STARID', cstarid)
    if cstarid gt 0 then begin
        objname = zspec.starid
        zspec = struct_addtags(replicate({objname:''}, n_elements(zspec)), zspec)
        zspec.objname = zspec.starid
        zspec = struct_trimtags(zspec, except_tags='STARID')
    endif
    wzspec1d = where(strtrim(tags, 2) eq 'ZSPEC1D_FILE', czs1d)
    if czs1d eq 0 then begin
        newstr = replicate({zspec1d_file:''}, n_elements(zspec))
        zspec = struct_addtags(zspec, newstr)
    endif


    ;read the zspec notes file
    newstr = replicate({zspec_slitid:-1, zspec_comment:' '}, n_elements(zspec))
    zspec = struct_addtags(zspec, newstr)
    notesfile = '/home/pemla/raja/m31_deimos/data_status/redux/s08/SAVE/'+mask+'.notes.txt'
    if file_test(notesfile) then begin
        comment = ' '
        zspec_slitid = -1L
        objno = ' '

        openr, 1, notesfile
        while ~eof(1) do begin
            point_lun, -1, pos
            readf, 1, comment, format='(2X,A1)'
            if comment eq '0' or comment eq '1' or comment eq '2' or comment eq '3' or comment eq '4' or comment eq '5' or comment eq '6' or comment eq '7' or comment eq '8' or comment eq '9' then begin
                point_lun, 1, pos
                readf, 1, zspec_slitid, objno, comment, format='(I3,14X,A9,2X,A80)'
                if strmid(strtrim(objno, 2), 0, 3) eq 'ser' then continue
                w = where(long(zspec.objname) eq long(objno), c)
                if c gt 0 then begin
                    zspec[w].zspec_slitid = zspec_slitid
                    zspec[w].zspec_comment = comment
                endif
            endif
        endwhile
        close, 1
    endif

    ;prepend spec1dpath to relevant filenames
    ;print, strmid(zspec[0].spec1d_file, 0, 4)
    if strmid(zspec[0].spec1d_file, 1, 4) eq '/net' then begin
        zspec.spec1d_file = strtrim(zspec.spec1d_file, 2)
        zspec.zspec1d_file = czs1d gt 0 ? strtrim(zspec.zspec1d_file, 2) : 'none'
        zspec.spec2d_file = strtrim(zspec.spec2d_file, 2)
    endif else begin
        if strmatch(spec1dpath, '*NEW_SPEC1D*') then begin
            for i=0,n_elements(zspec)-1 do begin
                zspec[i].spec1d_file = strtrim(spec1dpath+(strsplit(zspec[i].spec1d_file, '/', /extract))[1], 2)
                zspec[i].zspec1d_file = czs1d gt 0 ? strtrim(spec1dpath+(strsplit(zspec[i].zspec1d_file, '/', /extract))[1], 2) : 'none'
            endfor
            zspec.spec2d_file = strtrim('/net/pemla/data1/data2/m31_deimos/work/'+zspec.spec2d_file, 2)
        endif else begin
            zspec.spec1d_file = strtrim(spec1dpath+zspec.spec1d_file, 2)
            zspec.zspec1d_file = czs1d gt 0 ? strtrim(spec1dpath+zspec.zspec1d_file, 2) : 'none'
            zspec.spec2d_file = strtrim(spec1dpath+zspec.spec2d_file, 2)
        endelse
    endelse
    if strmatch(spec1dpath, '*f05_mgeha*') then begin
        for i=0,n_elements(zspec)-1 do begin
            zspec[i].spec1d_file = (strsplit(zspec[i].spec1d_file, '.gz', /extract, /regex))[0]
            zspec[i].zspec1d_file = czs1d gt 0 ? (strsplit(zspec[i].zspec1d_file, '.gz', /extract, /regex))[0] : 'none'
        endfor
    endif
    ;print, 'zspecfile: '+zspecfile
    nzspec = n_elements(zspec)    

    newzstr = replicate({vraw:-999d, verr:-999d, vhelio:-999d, vcorr:-999d}, nzspec)
    zspec = struct_trimtags(zspec, except_tags=['VRAW', 'VERR', 'VHELIO', 'VCORR'])
    zspec = struct_addtags(zspec, newzstr)

    c = 2.99792458d5

    goodfile=bytarr(nzspec)+1

    for i=0,nzspec-1 do begin
        if ~file_test(zspec[i].spec1d_file) then begin 
            message, zspec[i].spec1d_file+' not found', /info
            goodfile[i]=0
            spec = mrdfits(strcompress(zspec[i].spec1d_file, /rem), 1, hdr, /silent)
            mjd = sxpar(hdr, 'MJD-OBS') 
            jd = double(mjd) + 2400000.5
        endif

        type = size(zspec[i].ra, /type)
        case type of
            5: begin
                ra = zspec[i].ra
                dec = zspec[i].dec
            end
            7: begin
                get_coords, coords, instring=zspec[i].ra+'  '+zspec[i].dec
                ra = coords[0]*15d
                dec = coords[1]
            end
            else: message, 'Unknown RA/Dec format'
        endcase
        zspec[i].vhelio = (helio_deimos(ra, dec, jd=jd))[0]
        zspec[i].vraw = c*zspec[i].z
        zspec[i].verr = c*zspec[i].z_err
        zspec[i].vcorr = zspec[i].vraw - c*zspec[i].aband - zspec[i].vhelio
        ;zspec[i].v0 = c*(zspec[i].z) - zspec[i].vhelio
    endfor
    


    ;create the Na EW structure template
    na = {na_ew:-999d, $
          na_ewerr:-999d}
    na = replicate(na, nzspec)

    ;calculate Na EW
    wgoodfile=where(goodfile eq 1, cgoodfile)
    for i=0,cgoodfile-1 do begin
        na[wgoodfile[i]].na_ew = na_ew(zspec[wgoodfile[i]].spec1d_file, zspec[wgoodfile[i]].z, err=naerr, stop=long(zspec[wgoodfile[i]].objname) eq 4002236 ? 1 : 0)
        na[wgoodfile[i]].na_ewerr = naerr
        ;print, string(27B) + '[1A' + string(double(i)/double(nzspec-1)*100., format='(D5.1)')+'%'
    endfor
    
    ;create the Ca EW structure template
    ca = {sig_ca:-999d, $
          sig_ca_err:-999d, $
;          ;ca_a_ew:0d, $
;          ;ca_b_ew:0d, $
;          ;ca_c_ew:0d, $
;          ;ca_a_ew_err:0d, $
;          ;ca_b_ew_err:0d, $
;          ;ca_c_ew_err:0d, $
;          ;cont:0d, $
;          ;cont_err:0d
          feh_spec:-999d $
          }
    ca = replicate(ca, nzspec)

    ;calculate Ca EW
    for i=0,cgoodfile-1 do begin
        ca[wgoodfile[i]].sig_ca = ca_ew(zspec[wgoodfile[i]].spec1d_file, zspec[wgoodfile[i]].z, err=caerr);, stop=long(zspec[wgoodfile[i]].objname) eq 4002236 ? 1 : 0)
        ca[wgoodfile[i]].sig_ca_err = caerr
        ;print, string(27B) + '[1A' + string(double(i)/double(nzspec-1)*100., format='(D5.1)')+'%'
    endfor

    ;read the ASCII photometry file
    case field of
        'a19': photfileformat = 'L,X,X,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,X,D,D,D,D,X'
        'b15': photfileformat = 'L,X,X,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,X,D,D,D,D,X'
        'm11': photfileformat = 'L,X,X,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,X,D,D,D,D,X'
        else: photfileformat = 'L,X,X,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,X,D,D,D,D'
    endcase
    case photformat of
        'prob': readcol, photfile, id, D, Derr, M, Merr, T, Terr, chi, round, flag, prob, l, b, EBV, AM, EMT, EMD, xsi, eta, radii, gprob, gerr, ierr, serr, format=photfileformat, /silent, comment='#'
        'gprob': begin
            readcol, photfile, id, D, Derr, M, Merr, T, Terr, chi, round, flag, prob, radii, l, b, EBV, EMT, EMD, AM, gprob, gerr, ierr, serr, format='L,X,X,X,X,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent, comment='#'
            xsi = replicate(-999, n_elements(id))
            eta = replicate(-999, n_elements(id))
        end
        'phot': begin
            readcol, photfile, id, D, Derr, M, Merr, T, Terr, chi, round, l, b, EBV, EMT, EMD, AM, format='L,X,X,X,X,D,D,D,D,D,D,D,D,X,D,D,D,D,D,D', /silent, comment='#' ;round is actually sharp!!
            flag = replicate(-999, n_elements(id))
            prob = replicate(-999d, n_elements(id))
            xsi = replicate(-999d, n_elements(id))
            eta = replicate(-999d, n_elements(id))
            radii = replicate(-999d, n_elements(id))
            gprob = replicate(-999d, n_elements(id))
            gerr = replicate(-999d, n_elements(id))
            ierr = replicate(-999d, n_elements(id))
            serr = replicate(-999d, n_elements(id))
        end
        'cat': begin
            readcol, photfile, id, ra, dec, d, derr, t, terr, m, merr, chi, round, format='L,D,D,D,D,D,D,D,D,D,D', /silent, comment='#' ;round is actually sharp!!
            flag = replicate(-999, n_elements(id))
            prob = replicate(-999d, n_elements(id))
            l = replicate(-999d, n_elements(id))
            b = replicate(-999d, n_elements(id))
            ebv = replicate(-999d, n_elements(id))
            am = replicate(-999d, n_elements(id))
            emt = replicate(-999d, n_elements(id))
            emd = replicate(-999d, n_elements(id))
            xsi = replicate(-999d, n_elements(id))
            eta = replicate(-999d, n_elements(id))
            radii = replicate(-999d, n_elements(id))
            gprob = replicate(-999d, n_elements(id))
            gerr = replicate(-999d, n_elements(id))
            ierr = replicate(-999d, n_elements(id))
            serr = replicate(-999d, n_elements(id))
        end
    endcase

    ;create the photometry structure template
    phot = {idphot:-999L, $
            dm:-999d, $
            d:-999d, $
            derr:-999d, $
            m:-999d, $
            merr:-999d, $
            t:-999d, $
            terr:-999d, $
            i0:-999d, $
            i0err:-999d, $
            vi0:-999d, $
            vi0err:-999d, $
            xcmd:-999d, $
            xcmderr:-999d, $
            ycmd:-999d, $
            ycmderr:-999d, $
            feh_phot:-999d, $
            feh_photerr:-999d, $
            chi:-999d, $
            round:-999d, $
            l:-999d, $
            b:-999d, $
            ebv:-999d, $
            am:-999d, $
            emt:-999d, $
            emd:-999d, $
            xsi:-999d, $
            eta:-999d, $
            radii:-999d, $
            gprob:-999d, $
            flag_uva:-999d, $
            morphprob_uva:-999d, $
            gerr_uva:-999d, $
            ierr_uva:-999d, $
            serr_uva:-999d}
    phot = replicate(phot, n_elements(id))

    case savflag of
        'n147': dm = 24.10
        'n185': dm = 23.89
        'n205': dm = 24.58
        'and1': dm = 24.3608
        'and2': dm = 24.0712
        'and3': dm = 24.3724
        'and10': dm = 24.33
        'and14': dm = 24.7
        'm31': dm = 24.4688
    endcase

    ;fill the values
    phot.idphot = id
    phot.dm = dm
    phot.d = d
    phot.derr = derr
    phot.m = m
    phot.merr = merr
    phot.t = t
    phot.terr = terr
    phot.chi = chi
    phot.round = round
    phot.flag_uva = flag
    phot.morphprob_uva = prob
    phot.l = l
    phot.b = b
    phot.ebv = ebv
    phot.am = am
    phot.emt = emt
    phot.emd = emd
    phot.xsi = xsi
    phot.eta = eta
    phot.radii = radii
    phot.gprob = gprob
    phot.gerr_uva = gerr
    phot.ierr_uva = ierr
    phot.serr_uva = serr

    ;convert M,T to V,I
    if photfileformat eq 'cat' then begin
        phot.i0 = t
        phot.i0err = terr
        phot.vi0 = -0.006 + 0.8*(m-t)
        phot.vi0err = 0.8*sqrt(merr^2. + terr^2.)
    endif else begin
        phot.i0 = t - am + emt
        phot.i0err = terr
        phot.vi0 = -0.006 + 0.8*((m-t)-emt)
        phot.vi0err = 0.8*sqrt(merr^2. + terr^2.)
    endelse

    ;calculate X_CMD, Y_CMD, and [Fe/H]_phot
    xyfeh_arg, phot.i0, phot.vi0, feh, x, y, phot.i0err, phot.vi0err, feherr, xerr, yerr, /v_i, /van, savflag=savflag
    phot.feh_phot = feh
    phot.xcmd = x
    phot.ycmd = y
    phot.feh_photerr = feherr
    phot.xcmderr = xerr
    phot.ycmderr = yerr
    

    ;match zspec and phot structures by object ID
    match, long(zspec.objname), phot.idphot, w1, w2, count=c
    print, 'Found '+strtrim(c, 2)+' of '+strtrim(nzspec, 2)+' zspec objects in the photometry catalog.'
    s = sort(w1)
    zspec = zspec[w1[s]]
    na = na[w1[s]]
    ca = ca[w1[s]]
    phot = phot[w2[s]]
    if photfileformat eq 'cat' then begin
        euler, zspec.ra, zspec.dec, l, b, 1
        phot.l = l
        phot.b = b
        ebv = dust_getval(phot.l, phot.b, /interp)
        phot.vi0 -= 3.1*ebv - 0.482*3.1*ebv
        phot.i0 -= 0.482*3.1*ebv        
    endif

    ;concatenate all three structures
    alldata = struct_addtags(zspec, na)
    alldata = struct_addtags(alldata, ca)
    alldata = struct_addtags(alldata, phot)

    ;calculate [Fe/H]_spec
    v_HB = dm + 0.88                   ;Clem, J. L. 2005, Ph.D. thesis, Univ. Victoria
    feh_spec = -2.66 + 0.42 * (alldata.sig_ca - 0.64 * (v_HB - (alldata.vi0 + alldata.i0)))
    alldata.feh_spec = feh_spec
    ;print, feh_spec

    ;get rid of unwanted tags
    alldata = struct_trimtags(alldata, except_tags=['IDPHOT', 'IMAG', 'RMAG'])          ;IMAG and RMAG are 0's anyway

    ;calculate likelihood and class
    newstr = replicate({likelihood:-999d, m31class:0}, n_elements(alldata))
    alldata = struct_addtags(alldata, newstr)
    
    alldata.likelihood = m31prob(alldata.vcorr, alldata.gprob, alldata.vi0, alldata.na_ew, alldata.xcmd, alldata.ycmd, alldata.feh_phot, alldata.feh_spec, class=class, savflag=savflag)
    alldata.m31class = class

    ;write the alldata structure to a binary FITS file
    alldatafile = alldatapath+mask+'.fits'
    mwrfits, alldata, alldatafile, /create

    ;gzip the alldata FITS file
    spawn, 'gzip -f '+alldatafile
end
