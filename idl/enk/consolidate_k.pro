pro consolidate, mask, zspecpath=zspecpath, spec1dpath=spec1dpath, alldatapath=alldatapath
;consolidate, 'a19_1', alldatapath='/net/pemla/data1/data2/m31_deimos/contest/'
    ;default zspecpath
    if ~keyword_set(zspecpath) then zspecpath = '/net/pemla/data1/data2/m31_deimos/work/zresult/'
    if strmid(zspecpath, 0, 1, /reverse_offset) ne '/' then zspecpath += '/'

    ;default spec1dpath
    if ~keyword_set(spec1dpath) then spec1dpath = '/net/pemla/data1/data2/m31_deimos/work'
    if strmid(spec1dpath, 0, 1, /reverse_offset) eq '/' then spec1dpath = strmid(spec1dpath, 0, strlen(spec1dpath)-1)

    ;default alldatapath
    if ~keyword_set(alldatapath) then alldatapath = '/net/aegis2/u/priya/'
    ;if ~keyword_set(alldatapath) then alldatapath = '/net/schechter/a/ekirby/DEIMOS/phot/'
    ;if ~keyword_set(alldatapath) then alldatapath = '/net/pemla/data1/data3/kgilbert/DEIMOS/phot/'
    if strmid(alldatapath, 0, 1, /reverse_offset) ne '/' then alldatapath += '/'

    ;read the zspec FITS file
    zspecfile = zspecpath+'zspec.'+mask+'.fits'
    zspec = mrdfits(zspecfile, 1)

    ;prepend spec1dpath to relevant filenames
    ;print, strmid(zspec[0].spec1d_file, 0, 4)
    if strmid(zspec[0].spec1d_file, 1, 4) eq '/net' then begin
        zspec.spec1d_file = strtrim(zspec.spec1d_file, 2)
        zspec.zspec1d_file = strtrim(zspec.zspec1d_file, 2)
        zspec.spec2d_file = strtrim(zspec.spec2d_file, 2)
    endif else begin
        zspec.spec1d_file = strtrim(spec1dpath+zspec.spec1d_file, 2)
        zspec.zspec1d_file = strtrim(spec1dpath+zspec.zspec1d_file, 2)
        zspec.spec2d_file = strtrim(spec1dpath+zspec.spec2d_file, 2)
    endelse
    print, 'zspecfile: '+zspecfile
    nzspec = n_elements(zspec)    

    ;create the Na EW structure template
    na = {na_ew:0d, $
          na_ewerr:0d}
    na = replicate(na, nzspec)

    ;calculate Na EW
    for i=0,nzspec-1 do begin
        na[i].na_ew = na_ew(zspec[i].spec1d_file, zspec[i].z, err=naerr, stop=long(zspec[i].objname) eq 4002236 ? 1 : 0)
        na[i].na_ewerr = naerr
        ;print, string(27B) + '[1A' + string(double(i)/double(nzspec-1)*100., format='(D5.1)')+'%'
    endfor
    
    ;create the Ca EW structure template
    ca = {sig_ca:0d, $
          sig_ca_err:0d, $
;          ;ca_a_ew:0d, $
;          ;ca_b_ew:0d, $
;          ;ca_c_ew:0d, $
;          ;ca_a_ew_err:0d, $
;          ;ca_b_ew_err:0d, $
;          ;ca_c_ew_err:0d, $
;          ;cont:0d, $
;          ;cont_err:0d
          feh_spec:0d }
    ca = replicate(ca, nzspec)

    ;calculate Ca EW
    for i=0,nzspec-1 do begin
        ca[i].sig_ca = ca_ew(zspec[i].spec1d_file, zspec[i].z, err=caerr, stop=long(zspec[i].objname) eq 4002236 ? 1 : 0)
        ca[i].sig_ca_err = caerr
        ;print, string(27B) + '[1A' + string(double(i)/double(nzspec-1)*100., format='(D5.1)')+'%'
        print, 'ca', i
    endfor

    ;extract the field name from the mask name
    field = (strsplit(mask, '_', /extract))[0]

    ;read the ASCII photometry file
    ;photfile = '/net/pemla/u/raja/m31_deimos/data_status/redux_f07/work/phot/all.'+field+'.prob'
    ;photfileformat = field eq 'm8' ? 'L,X,X,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,X,D,D,D,D' : 'L,X,X,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,X,D,D,D,D,X'
    ;readcol, photfile, id, D, Derr, M, Merr, T, Terr, chi, round, flag, prob, l, b, EBV, AM, EMT, EMD, xsi, eta, radii, gprob, gerr, ierr, serr, format=photfileformat, /silent, comment='#'
    photfile = '/net/pemla/data/home_overflow/HSTCY15/m31_etc_acs/ppr2/clean_sample/deimosmatched2.dat_crop'
    ;photfileformat = 'L,L,D,D,X,X,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,L,L,D,D,D,D,D,D,D,L'
    photfileformat = 'L,L,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,L,L'
    readcol, photfile, count, photid, r, dec, D, Derr, M, Merr, T, Terr, chi, round, l, b, EBV, EMT, EMD, AM, gprob, gerr, g3, g4, id, slit, format=photfileformat, /silent, comment='#'
    flag=9999
    prob=9999
    xsi=9999
    eta=9999
    radii=9999
    ierr=9999
    serr=9999

    ;create the photometry structure template
    phot = {idphot:-999L, $
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

    ;fill the values
    phot.idphot = id
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
    phot.i0 = t - am + emt
    phot.i0err = terr
    phot.vi0 = -0.006 + 0.8*((m-t)-emt)
    phot.vi0err = 0.8*sqrt(merr^2. + terr^2.)

    ;calculate X_CMD, Y_CMD, and [Fe/H]_phot
    xyfeh_arg, phot.i0, phot.vi0, feh, x, y, phot.i0err, phot.vi0err, feherr, xerr, yerr, /v_i, /van, savflag='m31'
    phot.feh_phot = feh
    phot.xcmd = x
    phot.ycmd = y
    phot.feh_photerr = feherr
    phot.xcmderr = xerr
    phot.ycmderr = yerr
    
    print, 'done with xyfeh'
    
    ;print, feh_spec

    ;match zspec and phot structures by object ID
    match, long(zspec.objname), phot.idphot, w1, w2, count=c
    print, 'Found '+strtrim(c, 2)+' of '+strtrim(nzspec, 2)+' zspec objects in the photometry catalog.'
    s = sort(w1)
    zspec = zspec[w1[s]]
    na = na[w1[s]]
    ca = ca[w1[s]]
    phot = phot[w2[s]]

   
    print, 'now concatenating'

    ;concatenate all three structures
    alldata = struct_addtags(zspec, na)
    alldata = struct_addtags(alldata, ca)
    alldata = struct_addtags(alldata, phot)

    ;calculate [Fe/H]_spec
    v_HB = 25.17
    f_spec = -2.66 + 0.42 * (alldata.sig_ca - 0.64 * (v_HB - (alldata.vi0 + alldata.i0)))
    alldata.feh_spec = f_spec
 
    ;get rid of unwanted tags
    alldata = struct_trimtags(alldata, except_tags=['IDPHOT', 'IMAG', 'RMAG'])          ;IMAG and RMAG are 0's anyway

    ;write the alldata structure to a binary FITS file
    alldatafile = alldatapath+mask+'.fits'
    mwrfits, alldata, alldatafile, /create

    print, 'gzipping'
    ;gzip the alldata FITS file
    spawn, 'gzip '+alldatafile
end
