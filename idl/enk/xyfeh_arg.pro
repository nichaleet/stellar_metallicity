pro xyfeh_arg, mag_in, clr_in, feh_out, x_out, y_out, err_mag_in, err_clr_in, err_feh, err_x, err_y, r_i=r_i, v_i=v_i, van=van, padova=padova, savflag=savflag

    if keyword_set(van)+keyword_set(padova) ne 1 then message, 'You must choose VandenBerg or Padova isochrones.'
    ;if ~keyword_set(infile) then message, 'You must specify an input file.'
    if keyword_set(v_i)+keyword_set(r_i) ne 1 then message, 'You must specify one and only one color (V-I or R-I).'

    if ~keyword_set(savflag) then savflag = ''
    savflag = '_'+savflag
    if keyword_set(v_i) then savflag = savflag+'_vi'
    if keyword_set(r_i) then savflag = savflag+'_ri'
    if keyword_set(van) then savflag = savflag+'_van'
    if keyword_set(padova) then savflag = savflag+'_padova'
    restore, getenv('M31')+'photmetal/xyfeh'+savflag+'.sav'

    n = n_elements(mag_in)

    feh_out = interp2d(feh_grid, clr_grid, mag_grid, clr_in, mag_in, missing=-999., /regular)
    x_out = interp2d(x_grid, clr_grid, mag_grid, clr_in, mag_in, missing=-999., /regular)
    y_out = interp2d(y_grid, clr_grid, mag_grid, clr_in, mag_in, missing=-999., /regular)
    
    if n_params() gt 5 then begin
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
end
