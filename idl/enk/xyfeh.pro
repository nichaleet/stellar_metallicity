pro xyfeh, infile=infile, outfile=outfile, r_i=r_i, v_i=v_i, van=van, padova=padova, master=master, ricngood=ricngood, ricnradec=ricnradec, karrie=karrie, merge=merge, cut2=cut2, out2=out2, m31_phot=m31_phot, besancon=besancon, plot=plot, savflag=savflag, jason070509=jason070509, de=de, dsph=dsph, alldata=alldata

    if keyword_set(van)+keyword_set(padova) ne 1 then message, 'You must choose VandenBerg or Padova isochrones.'
    ;if ~keyword_set(infile) then message, 'You must specify an input file.'
    if keyword_set(v_i)+keyword_set(r_i) ne 1 then message, 'You must specify one and only one color (V-I or R-I).'
    if keyword_set(master)+keyword_set(ricngood)+keyword_set(ricnradec)+keyword_set(karrie)+keyword_set(merge)+keyword_set(cut2)+keyword_set(out2)+keyword_set(m31_phot)+keyword_set(besancon)+keyword_set(jason070509)+keyword_set(de)+keyword_set(dsph)+keyword_set(alldata) ne 1 then message, 'You must specify one and only one file format keyword.'

    if keyword_set(infile) then infile = keyword_set(m31_phot) or keyword_set(besancon) or keyword_set(alldata) ? infile : '/net/schechter/a/ekirby/m31_sandbox/'+infile else infile = dialog_pickfile(path='/net/schechter/a/ekirby/m31_sandbox/', title='Please select an input photometry file')
    if ~keyword_set(outfile) then outfile = '/net/schechter/a/ekirby/m31_sandbox/'+file_basename(infile)+'.xyfeh'

    if ~keyword_set(savflag) then savflag = ''
    savflag = '_'+savflag
    if keyword_set(v_i) then savflag = savflag+'_vi'
    if keyword_set(r_i) then savflag = savflag+'_ri'
    if keyword_set(van) then savflag = savflag+'_van'
    if keyword_set(padova) then savflag = savflag+'_padova'
    restore, getenv('M31')+'photmetal/xyfeh'+savflag+'.sav'


    ;format for MASTER files: I, V-I
    if keyword_set(master) then readcol, infile, in_1, in_2, field_in, in_4, in_5, mag_in, clr_in, in_8, x_in, y_in, feh_phot_in, feh_spec_in, in_13, li_in, format='I,L,A,F,F,F,F,F,F,F,F,F,F,F'

    ;format for RICNgood files: I, R-I
    if keyword_set(ricngood) then begin
        readcol, infile, id, ra_h, ra_m, ra_s, dec_d, dec_m, dec_s, mag_in, err_mag_in, clr_in, err_clr_in, cntio, err_cntio, format='L,F,F,F,F,F,F,F,F,F,F,F,F', skipline=3
        ra = (ra_h + ra_m/60. + ra_s/3600.)*360./24.
        dec = dec_d + dec_m/60. + dec_s/3600.
    endif

    ;format for RICNradec files: I, R-I
    if keyword_set(ricnradec) then begin
        readcol, infile, id, ra, dec, mag_in, err_mag_in, clr_in, err_clr_in, cntio, err_cntio, format='L,A,A,F,F,F,F,F,F'
        ra = '0'+ra
        dec = '+'+dec
        coords = sex2deg([[ra], [dec]])
        ra = coords[*,0]
        dec = coords[*,1]
    endif

    ;format for Karrie's abridged MASTER file: I, V-I
    if keyword_set(karrie) then readcol, infile, id, mask, vel, fdd0, mag_in, clr_in, ew_na, feh_spec, format='L,A,F,F,F,F,F,F'
    
    ;format for merge files
    ;if keyword_set(merge) then readcol, infile, field_in, id_in, i_in, rasex_in, decsex_in, ra_in, dec_in, in_8, in_9, in_10, in_11, in_12, in_13, mag_in, clr_in, format='A,L,L,A,A,F,F,F,F,F,F,F,L,F,F'
     if keyword_set(merge) then readcol, infile, field_in, id_in, i_in, rasex_in, decsex_in, ra_in, dec_in, in_8, in_9, in_10, in_11, in_12, in_13, in_14, mag_in, err_mag_in, clr_in, err_clr_in, format='A,L,L,A,A,F,F,F,F,F,F,F,L,F,F,F,F,F'

    ;format for cut2 files: I, V-I
    if keyword_set(cut2) then readcol, infile, field_in, id_in, i_in, rasex_in, decsex_in, ra_in, dec_in, mag_in, clr_in, in_10, in_11, in_12, in_13, li_in, in_15, format='A,L,L,A,A,F,F,F,F,F,F,F,F,L,L'

    ;format for out2 files: I, V-I
    if keyword_set(out2) then readcol, infile, field_in, id_in, in_3, rasex_in, decsex_in, ra_in, dec_in, mag_in, clr_in, in_10, in_11, in_12, in_13, in_14, li_in, i_in, format='A,L,L,A,A,F,F,F,F,F,F,F,F,F,L,L'

    ;format for m31_phot.fits: I, V-I
    if keyword_set(m31_phot) then begin
        phot = mrdfits(infile, 1)
        mag_in = phot.t
        clr_in = phot.v - phot.t
        if m31_phot eq 2 or m31_phot eq 4 then clr_in = clr_in - (phot.merr > phot.terr)
        if m31_phot eq 3 or m31_phot eq 5 then clr_in = clr_in + (phot.merr > phot.terr)
    endif

    if keyword_set(besancon) then begin
        besancon = mrdfits(infile, 1)
        mag_in = besancon.t
        clr_in = besancon.vi
    endif
        
    ;format for jason.070509: I, V-I
    if keyword_set(jason070509) then readcol, infile, mag_in, clr_in, format='D,D'

    ;format for dEs: I, R-I
    if keyword_set(de) then readcol, infile, id_in, ra_in, dec_in, imagpriority_in, junk1_in, junk2_in, junk3_in, junk4_in, mag_in, err_mag_in, clr_in, vi_0_in, err_clr_in, mask_in, radeg_in, decdeg_in, vel_in, velcorr_in, velerr_in, S2N_in, qual_in, format='L,A,A,D,I,I,I,I,D,D,D,D,D,A,D,D,D,D,D,D,I'

    ;format for dSphs: I, V-I
    if keyword_set(dsph) then readcol, infile, id_in, ra_in, dec_in, imagpriority_in, junk1_in, junk2_in, junk3_in, mag_in, err_mag_in, clr_in, err_clr_in, mask_in, radeg_in, decdeg_in, vel_in, velcorr_in, velerr_in, S2N_in, qual_in, format='L,A,A,D,I,I,I,D,D,D,D,A,D,D,D,D,D,D,I'

    ;format for alldata: I, V-I
    if keyword_set(alldata) then begin
        readcol, infile, maskname, id, slit, ras, decs, ra, dec, vraw, vhelio, vcorr, verr, sn, zq, d, derr, m, merr, t, terr, chi, round, flag, prob, l, b, $
                 format='A,L,L,A,A,D,D,D,D,D,D,D,I,X,X,X,D,D,D,D,D,D,D,D,D,D,D,D,X,X,X,X,X,X,X,X,X,X,X,X,X'
        readcol, infile, ebv, am, emt, emd, xsi, eta, radii, gprob, gerr, ierr, serr, format='X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,D,D,D,D,D,D,D,X,D,D,D,D,X'
        m0 = m - am
        mt0 = (m-t)-emt
        err_mt0 = sqrt((merr)^2. + (terr)^2.)
        mag_in = m0 - mt0
        clr_in = -0.006 + 0.8*mt0
        err_mag_in = terr
        err_clr_in = 0.8*abs(err_mt0)
    endif

    n = n_elements(mag_in)

    feh_out = interp2d(feh_grid, clr_grid, mag_grid, clr_in, mag_in, missing=-999., /regular)
    x_out = interp2d(x_grid, clr_grid, mag_grid, clr_in, mag_in, missing=-999., /regular)
    y_out = interp2d(y_grid, clr_grid, mag_grid, clr_in, mag_in, missing=-999., /regular)
    

    if keyword_set(m31_phot) then begin
        phot.x[m31_phot-1] = x_out
        phot.y[m31_phot-1] = y_out
        phot.feh[m31_phot-1] = feh_out
        mwrfits, phot, infile, /create
    endif
    
    if keyword_set(besancon) then begin
        besancon.x = x_out
        besancon.y = y_out
        besancon.feh = feh_out
        mwrfits, besancon, infile, /create
    endif
    
    if keyword_set(ricnradec) or keyword_set(ricngood) or keyword_set(de) or keyword_set(dsph) or keyword_set(merge) or keyword_set(alldata) then begin
       remember = 0
       nmc = 1000
       err_feh = fltarr(n)
       err_x = fltarr(n)
       err_y = fltarr(n)
       if remember eq 1 then begin
          feh_mc = fltarr(nmc, n)
          x_mc = fltarr(nmc, n)
          y_mc = fltarr(nmc, n)
       endif
       print, 'interpolation and error calculation progress:'
       print, ' '
       for i=0L,n-1 do begin
          seed = float(strmid(string(systime(/seconds), format='(D11.0)'), 5, 5, /reverse_offset))
          rand = randomn(seed+i, nmc, 2)
          clr_mc = clr_in[i]+rand[*,0]*err_clr_in[i]
          mag_mc = mag_in[i]+rand[*,1]*err_mag_in[i]
          if remember eq 1 then begin
             feh_mc[*,i] = interp2d(feh_grid, clr_grid, mag_grid, clr_mc, mag_mc, missing=-999., /regular)
             x_mc[*,i] = interp2d(x_grid, clr_grid, mag_grid, clr_mc, mag_mc, missing=-999., /regular)
             y_mc[*,i] = interp2d(y_grid, clr_grid, mag_grid, clr_mc, mag_mc, missing=-999., /regular)
             err_feh[i] = stddev(feh_mc[*,i])
             err_x[i] = stddev(x_mc[*,i])
             err_y[i] = stddev(y_mc[*,i])
          endif else begin
             feh_mc = interp2d(feh_grid, clr_grid, mag_grid, clr_mc, mag_mc, missing=-999., /regular)
             x_mc = interp2d(x_grid, clr_grid, mag_grid, clr_mc, mag_mc, missing=-999., /regular)
             y_mc = interp2d(y_grid, clr_grid, mag_grid, clr_mc, mag_mc, missing=-999., /regular)
             err_feh[i] = stddev(feh_mc)
             err_x[i] = stddev(x_mc)
             err_y[i] = stddev(y_mc)             
          endelse
          print, string(27B) + '[1A' + string(double(100*(i+1)) / double(n), format='(D5.1)') + '%'
       endfor
    endif

    print, ' '
    print, 'writing output to file ...'
    if keyword_set(clrtweak) then begin
        phot.feh = feh_out
    endif

    openw, 1, outfile, width='10000'
    for j=0L,n-1 do begin
        if keyword_set(master) then printf, 1, in_1[j], in_2[j], '   '+field_in[j], in_4[j], in_5[j], mag_in[j], clr_in[j], in_8[j], x_out[j], y_out[j], feh_out[j], feh_spec_in[j], in_13[j], li_in[j]
        if keyword_set(ricngood) or keyword_set(ricnradec) then printf, 1, id[j], ra[j], dec[j], mag_in[j], err_mag_in[j], clr_in[j], err_clr_in[j], cntio[j], err_cntio[j], feh_out[j], err_feh[j], x_out[j], err_x[j], y_out[j], err_y[j]
        if keyword_set(karrie) then printf, 1, id[j], '   '+mask[j], vel[j], fdd0[j], mag_in[j], clr_in[j], ew_na[j], feh_spec[j], feh_out[j], x_out[j], y_out[j]
        if keyword_set(merge) then printf, 1, '   '+field_in[j], id_in[j], i_in[j], '   '+rasex_in[j], '   '+decsex_in[j], ra_in[j], dec_in[j], in_8[j], in_9[j], in_10[j], in_11[j], in_12[j], in_13[j], in_14[j], mag_in[j], err_mag_in[j], clr_in[j], err_clr_in[j], feh_out[j], err_feh[j], x_out[j], err_x[j], y_out[j], err_y[j]
        if keyword_set(cut2) then printf, 1, '   '+field_in[j], id_in[j], i_in[j], '   '+rasex_in[j], '   '+decsex_in[j], ra_in[j], dec_in[j], mag_in[j], clr_in[j], in_10[j], in_11[j], in_12[j], in_13[j], li_in[j], in_15[j], feh_out[j], x_out[j], y_out[j]
        if keyword_set(out2) then printf, 1, '   '+field_in[j], id_in[j], in_3[j], '   '+rasex_in[j], '   '+decsex_in[j], ra_in[j], dec_in[j], mag_in[j], clr_in[j], in_10[j], in_11[j], in_12[j], in_13[j], in_14[j], li_in[j], i_in[j], feh_out[j], x_out[j], y_out[j]
        if keyword_set(jason070509) then printf, 1, mag_in[j], clr_in[j], feh_out[j], x_out[j], y_out[j]
        if keyword_set(de) then printf, 1, id_in[j], ra_in[j], dec_in[j], imagpriority_in[j], junk1_in[j], junk2_in[j], junk3_in[j], junk4_in[j], mag_in[j], err_mag_in[j], clr_in[j], vi_0_in[j], err_clr_in[j], mask_in[j], radeg_in[j], decdeg_in[j], vel_in[j], velcorr_in[j], velerr_in[j], S2N_in[j], qual_in[j], feh_out[j], err_feh[j], x_out[j], err_x[j], y_out[j], err_y[j]
        if keyword_set(dsph) then printf, 1, id_in[j], ra_in[j], dec_in[j], imagpriority_in[j], junk1_in[j], junk2_in[j], junk3_in[j], mag_in[j], err_mag_in[j], clr_in[j], err_clr_in[j], mask_in[j], radeg_in[j], decdeg_in[j], vel_in[j], velcorr_in[j], velerr_in[j], S2N_in[j], qual_in[j], feh_out[j], err_feh[j], x_out[j], err_x[j], y_out[j], err_y[j]
        if keyword_set(alldata) then printf, 1, maskname[j], id[j], slit[j], ras[j], decs[j], ra[j], dec[j], vraw[j], vhelio[j], vcorr[j], verr[j], sn[j], zq[j], d[j], derr[j], m[j], merr[j], t[j], terr[j], chi[j], round[j], flag[j], prob[j], l[j], b[j], ebv[j], am[j], emt[j], emd[j], xsi[j], eta[j], radii[j], gprob[j], gerr[j], ierr[j], serr[j], feh_out[j], err_feh[j], x_out[j], err_x[j], y_out[j], err_y[j]
    endfor
    close, 1

    if 0 then begin
        feh_phot = {field:' ', mag:-999., clr:-999., x_old:-999., x_new:-999., y_old:-999., y_new:-999., feh_spec:-999., feh_phot_old:-999., feh_phot_new:-999., li:-999.}
        feh_phot = replicate(feh_phot, n)
        feh_phot.field = field_in
        feh_phot.mag = mag_in
        feh_phot.clr = clr_in
        feh_phot.x_old = x_in
        feh_phot.x_new = x_out
        feh_phot.y_old = y_in
        feh_phot.y_new = y_out
        feh_phot.feh_spec = feh_spec_in
        feh_phot.feh_spec = feh_spec_in
        feh_phot.feh_phot_old = feh_phot_in
        feh_phot.feh_phot_new = feh_out
        feh_phot.li = li_in

        ;feh_phot = feh_phot[where(feh_phot.mag lt 50 and feh_phot.feh_phot_old gt -10 and feh_phot.feh_phot_old lt 5 and feh_phot.li gt 1 or (feh_phot.li eq 0 and in_13 gt 0.5))]
    endif

    if keyword_set(plot) then begin
        setplot
        !X.STYLE = 1
        !Y.STYLE = 1        

        w = where(feh_phot.x_old gt -5 and feh_phot.x_new gt -5 and feh_phot.x_old lt 5 and feh_phot.x_new lt 5)
        device, filename='x_compare.eps', /inches, xsize=6, ysize=6, xoffset=0.5, yoffset=0.5, /color, /encapsulated
        plot, feh_phot[w].x_old, feh_phot[w].x_new - feh_phot[w].x_old, xtitle='x!DJ!N', ytitle='x!DE!N - x!DJ!N', psym=1
        device, /close

        w = where(feh_phot.y_old gt -5 and feh_phot.y_new gt -5 and feh_phot.y_old lt 5 and feh_phot.y_new lt 5)
        device, filename='y_compare.eps', /inches, xsize=6, ysize=6, xoffset=0.5, yoffset=0.5, /color, /encapsulated
        plot, feh_phot[w].y_old, feh_phot[w].y_new - feh_phot[w].y_old, xtitle='y!DJ!N', ytitle='y!DE!N - y!DJ!N', psym=1
        device, /close

        w = where(feh_phot.feh_phot_old gt -10 and feh_phot.feh_phot_new gt -10 and feh_phot.feh_phot_old lt 5 and feh_phot.feh_phot_new lt 5)
        device, filename='feh_compare.eps', /inches, xsize=6, ysize=6, xoffset=0.5, yoffset=0.5, /color, /encapsulated
        plot, feh_phot[w].feh_phot_old, feh_phot[w].feh_phot_new - feh_phot[w].feh_phot_old, xtitle='[Fe/H]!DJ!N', ytitle='[Fe/H]!DE!N - [Fe/H]!DJ!N', psym=1, xrange=[-4, 0.7], yrange=[-0.5, 1.5]
        device, /close

        resetplot
    endif

    ;return, feh_phot
end
