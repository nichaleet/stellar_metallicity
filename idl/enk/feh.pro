;+
; NAME: feh
;
; PURPOSE: Calculate photometric metallicity with Victoria-Regina,
;    Yonsei-Yale, or Padova isochrones.
;
; CALLING SEQUENCE:
;    fehphot = feh(mag_in, clr_in, err_mag_in, err_clr_in, /color,
;       /isochrone, dm=dm, /alpha, /extrap)
;
; INPUTS:
;    mag_in:  An array of extinction-corrected magnitudes.
;    clr_in:  An array of dereddened colors.  mag_in is always the
;       redder filter. For example, if your color is V-I, then mag_in
;       must be I.  clr_in must have the same number of elements as
;       mag_in.
;    dm:      The absolute distance modulus (m-M)_0.  It should have
;       one element for all the input stars, or it should have as many
;       elements as the input stars.  In other words, each star can
;       have its own distance modulus.
;
; OPTIONAL INPUTS:
;    err_mag_in: An array of photometric errors on the magnitudes.  It
;       must have the same number of elements as mag_in.
;    err_clr_in: An array of photometric errors on the colors.  It
;       must have the same number of elements as mag_in.
;
; KEYWORD PARAMETERS:
;  Choose exactly one of the following color flags:
;    /ri:     Use R-I color and I magnitude.
;    /vi:     Use V-I color and I magnitude.
;    /bv:     Use B-V color and V magnitude.
;    /bi:     Use B-I color and I magnitude.
;    /br:     Use B-R color and R magnitude.
;    /vr:     Use V-R color and R magnitude.
;    /vj:     Use V-J color and J magnitude.
;    /vh:     Use V-H color and H magnitude.
;    /vk:     Use V-K color and K magnitude.
;    /jk:     Use J-K color and K magnitude.
;    /gisdss: Use SDSS g-i color and i magnitude.
;    /grsdss: Use SDSS g-r color and r magnitude.
;
;  Choose exactly one of the following isochrone flags:
;    /van:    Use Victoria-Regina (VandenBerg et al. 2006) isochrones.
;    /yy:     Use Yonsei-Yale (Demarque et al. 2004) isochrones.
;    /padova: Use Padova (Girardi et al. 2002) isochrones.
;
;  These are optional:
;    /alpha:  Use [alpha/Fe] = +0.3 instead of 0.0 isochrones.
;    /extrap: Extrapolate beyond the isochrone boundaries.
;
; OUTPUTS:
;  A structure with as many elements as mag_in.  The structure
;  contains the following tags:
;    AGES:    An array of stellar ages in Gyr.  The number of ages
;       depends on the choice of isochrone.
;    FEH:     An array of photometric metallicities corresponding to
;       each age.
;    TEFF:    An array of effective temperatures (K) corresponding to
;       each age.
;    LOGG:    An array of surface gravities (log cm/s^2) corresponding
;       to each age.
;    ERR_FEH: An array of errors on the photometric metallicity
;       corresponding to each element in FEH.  This tag exists only if
;       the inputs include photometric errors.
;    ERR_TEFF: An array of errors on the effective temperatures
;       corresponding to each element in TEFF.  This tag exists only
;       if the inputs include photometric errors.
;    ERR_LOGG: An array of errors on the surface gravities
;       corresponding to each element in LOGG.  This tag exists only
;       if the inputs include photometric errors.
;
; MODIFICATION HISTORY:
;    created Evan N. Kirby (UCSC, 2008 Aug 14)
;    ENK added B-R color, 2008 Sep 05
;    ENK added B-I color, 2012 Aug 08
;-


function feh, mag_in, clr_in, err_mag_in, err_clr_in, ri=ri, vi=vi, bv=bv, bi=bi, br=br, vr=vr, vj=vj, vh=vh, vk=vk, jk=jk, gisdss=gisdss, grsdss=grsdss, van=van, padova=padova, yy=yy, dm=dm, ddm=ddm, alpha=alpha, extrap=extrap
    dmin = dm

    if keyword_set(van)+keyword_set(padova)+keyword_set(yy) ne 1 then message, 'You must choose Victoria-Regina, Padova, or Yonsei-Yale isochrones.'
    if keyword_set(ri)+keyword_set(vi)+keyword_set(bv)+keyword_set(bi)+keyword_set(br)+keyword_set(vr)+keyword_set(vj)+keyword_set(vh)+keyword_set(vk)+keyword_set(jk)+keyword_set(gisdss)+keyword_set(grsdss) ne 1 then message, 'You must choose one and only one color.'
    if (keyword_set(vj)+keyword_set(vh)+keyword_set(vk)+keyword_set(jk) gt 0) and keyword_set(van) then message, 'Victoria-Regina isochrones do not support JHK.'

    noerror = n_params() eq 4 ? 0 : 1

    n = n_elements(mag_in)
    dm_uniq = dm[uniq(dm, sort(dm))]
    ndm = n_elements(dm_uniq)
    if n_elements(dm) eq 1 and n gt 1 then dm = replicate(dm, n)

    if keyword_set(ddm) and noerror eq 0 then err_mag_in = sqrt(err_mag_in^2. + ddm^2.)

    nmc = 1000
    if keyword_set(gisdss) then filestr = 'feh_gisdss_'
    if keyword_set(grsdss) then filestr = 'feh_grsdss_'
    if keyword_set(vi) then filestr = 'feh_vi_'
    if keyword_set(ri) then filestr = 'feh_ri_'
    if keyword_set(bv) then filestr = 'feh_bv_'
    if keyword_set(bi) then filestr = 'feh_bi_'
    if keyword_set(br) then filestr = 'feh_br_'
    if keyword_set(vr) then filestr = 'feh_vr_'
    if keyword_set(vj) then filestr = 'feh_vj_'
    if keyword_set(vh) then filestr = 'feh_vh_'
    if keyword_set(vk) then filestr = 'feh_vk_'
    if keyword_set(jk) then filestr = 'feh_jk_'
    alphastr = keyword_set(alpha) eq 1 ? 'a_' : '[0-9]'
    if keyword_set(van) then filestr = filestr+alphastr+'*_van.sav'
    if keyword_set(padova) then filestr = filestr+alphastr+'*_padova.sav'
    if keyword_set(yy) then filestr = filestr+alphastr+'*_yy.sav'
    savfiles = file_search(getenv('ahome')+'m31/photmetal/fehgrids', filestr)
    nsav = n_elements(savfiles)
    if noerror eq 0 then begin
        fehphot = {mag:-999., err_mag:-999., clr:-999., err_clr:-999., ages:fltarr(nsav)-999., feh:fltarr(nsav)-999., err_feh:fltarr(nsav)-999., teff:fltarr(nsav)-999., err_teff:fltarr(nsav)-999., logg:fltarr(nsav)-999., err_logg:fltarr(nsav)-999., logl:fltarr(nsav)-999., err_logl:fltarr(nsav)-999.}
    endif else begin
        fehphot = {mag:-999., clr:-999., ages:fltarr(nsav)-999., feh:fltarr(nsav)-999., teff:fltarr(nsav)-999., logg:fltarr(nsav)-999., logl:fltarr(nsav)-999.}
    endelse
    fehphot = replicate(fehphot, n)
    fehphot.mag = mag_in
    fehphot.clr = clr_in
    if noerror eq 0 then begin
        fehphot.err_mag = err_mag_in
        fehphot.err_clr = err_clr_in
        err_feh = fltarr(nsav, n)
        feh_out = fltarr(nsav, n)
    endif
    ages = fltarr(nsav)
    t = 0l
    ttot = n*nsav
    if noerror eq 0 then begin
        print, 'interpolation and error calculation progress:'
        print, ' '
    endif
    for ii=0,nsav-1 do begin
        restore, savfiles[ii]
        if keyword_set(extrap) then begin
            feh_grid = feh_grid_extrap
            logteff_grid = logteff_grid_extrap
            logg_grid = logg_grid_extrap
        endif
        fehphot.ages[ii] = replicate(age, n)

        for jj=0,ndm-1 do begin
            mag_grid += dm_uniq[jj]
            wdm = where(dm eq dm_uniq[jj], ndm_uniq)

            wgood = where(feh_grid ne -999 and logteff_grid ne -999 and logg_grid ne -999)

            feh_max = max(feh_grid[wgood])
            feh_min = min(feh_grid[wgood])
            feh_i = interp2d(feh_grid, clr_grid, mag_grid, clr_in[wdm], mag_in[wdm], missing=-999., /regular)
            fehphot[wdm].feh[ii] = feh_i
            w = where(feh_i le feh_min or feh_i ge feh_max, cw)
            if cw gt 0 then feh_i[w] = -999

            logteff_max = max(logteff_grid[wgood])
            logteff_min = min(logteff_grid[wgood])
            logteff_i = interp2d(logteff_grid, clr_grid, mag_grid, clr_in[wdm], mag_in[wdm], missing=-999., /regular)
            fehphot[wdm].teff[ii] = 10.0^(logteff_i)
            w = where(logteff_i le logteff_min or logteff_i ge logteff_max, cw)
            if cw gt 0 then logteff_i[w] = -999

            logg_max = max(logg_grid[wgood])
            logg_min = min(logg_grid[wgood])
            logg_i = interp2d(logg_grid, clr_grid, mag_grid, clr_in[wdm], mag_in[wdm], missing=-999., /regular)
            fehphot[wdm].logg[ii] = logg_i
            w = where(logg_i le logg_min or logg_i ge logg_max, cw)
            if cw gt 0 then logg_i[w] = -999

            log_l_max = max(log_l_grid[wgood])
            log_l_min = min(log_l_grid[wgood])
            log_l_i = interp2d(log_l_grid, clr_grid, mag_grid, clr_in[wdm], mag_in[wdm], missing=-999., /regular)
            fehphot[wdm].logl[ii] = log_l_i
            w = where(log_l_i le log_l_min or log_l_i ge log_l_max, cw)
            if cw gt 0 then log_l_i[w] = -999

            if noerror eq 0 then begin
                for i=0L,ndm_uniq-1 do begin
                    if feh_i[i] eq -999 then begin
                        fehphot[wdm[i]].err_feh[ii] = 0
                    endif else begin
                        seed = float(strmid(string(systime(/seconds), format='(D11.0)'), 5, 5, /reverse_offset))
                        rand = randomn(seed+i, nmc, 2)
                        clr_mc = clr_in[wdm[i]]+rand[*,0]*err_clr_in[wdm[i]]
                        mag_mc = mag_in[wdm[i]]+rand[*,1]*err_mag_in[wdm[i]]
                        if keyword_set(extrap) then begin
                            feh_mc = interp2d(feh_grid_extrap, clr_grid, mag_grid, clr_mc, mag_mc, /extrapolate, /regular)
                            wgood = where(feh_mc gt -5 and feh_mc lt 2.0, cgood)
                            if cgood ge 2 then fehphot[wdm[i]].err_feh[ii] = stddev(feh_mc[wgood]) else fehphot[wdm[i]].err_feh[ii] = -999
                        endif else begin
                            feh_mc = interp2d(feh_grid, clr_grid, mag_grid, clr_mc, mag_mc, /regular, missing=-999)
                            wgood = where(feh_mc gt -5 and feh_mc lt 2.0, cgood)
                            if cgood ge 2 then fehphot[wdm[i]].err_feh[ii] = stddev(feh_mc[wgood]) else fehphot[wdm[i]].err_feh[ii] = -999
                        endelse
                    endelse

                    if logteff_i[i] eq -999 then begin
                        fehphot[wdm[i]].err_teff[ii] = 0
                    endif else begin
                        seed = float(strmid(string(systime(/seconds), format='(D11.0)'), 5, 5, /reverse_offset))
                        rand = randomn(seed+i, nmc, 2)
                        clr_mc = clr_in[wdm[i]]+rand[*,0]*err_clr_in[wdm[i]]
                        mag_mc = mag_in[wdm[i]]+rand[*,1]*err_mag_in[wdm[i]]
                        if keyword_set(extrap) then begin
                            logteff_mc = interp2d(logteff_grid_extrap, clr_grid, mag_grid, clr_mc, mag_mc, /extrapolate, /regular)
                            wgood = where(10.^logteff_mc gt 3000 and 10.^logteff_mc lt 10000, cgood)
                            if cgood ge 2 then fehphot[wdm[i]].err_teff[ii] = abs(fehphot[wdm[i]].teff[ii]*alog(10.0))*stddev(logteff_mc[wgood]) else fehphot[wdm[i]].err_teff[ii] = -999
                        endif else begin
                            logteff_mc = interp2d(logteff_grid, clr_grid, mag_grid, clr_mc, mag_mc, /regular, missing=-999)
                            wgood = where(10.^logteff_mc gt 3000 and 10.^logteff_mc lt 10000, cgood)
                            if cgood ge 2 then fehphot[wdm[i]].err_teff[ii] = abs(fehphot[wdm[i]].teff[ii]*alog(10.0))*stddev(logteff_mc[wgood]) else fehphot[wdm[i]].err_teff[ii] = -999
                        endelse
                    endelse

                    if logg_i[i] eq -999 then begin
                        fehphot[wdm[i]].err_logg[ii] = 0
                    endif else begin
                        seed = float(strmid(string(systime(/seconds), format='(D11.0)'), 5, 5, /reverse_offset))
                        rand = randomn(seed+i, nmc, 2)
                        clr_mc = clr_in[wdm[i]]+rand[*,0]*err_clr_in[wdm[i]]
                        mag_mc = mag_in[wdm[i]]+rand[*,1]*err_mag_in[wdm[i]]
                        if keyword_set(extrap) then begin
                            logg_mc = interp2d(logg_grid_extrap, clr_grid, mag_grid, clr_mc, mag_mc, /extrapolate, /regular)
                            wgood = where(logg_mc gt -2.0 and logg_mc lt 7.0, cgood)
                            if cgood ge 2 then fehphot[wdm[i]].err_logg[ii] = stddev(logg_mc[wgood]) else fehphot[wdm[i]].err_logg[ii] = -999
                        endif else begin
                            logg_mc = interp2d(logg_grid, clr_grid, mag_grid, clr_mc, mag_mc, /regular, missing=-999)
                            wgood = where(logg_mc gt -2.0 and logg_mc lt 7.0, cgood)
                            if cgood ge 2 then fehphot[wdm[i]].err_logg[ii] = stddev(logg_mc[wgood]) else fehphot[wdm[i]].err_logg[ii] = -999
                        endelse
                    endelse

                    if log_l_i[i] eq -999 then begin
                        fehphot[wdm[i]].err_logl[ii] = 0
                    endif else begin
                        seed = float(strmid(string(systime(/seconds), format='(D11.0)'), 5, 5, /reverse_offset))
                        rand = randomn(seed+i, nmc, 2)
                        clr_mc = clr_in[wdm[i]]+rand[*,0]*err_clr_in[wdm[i]]
                        mag_mc = mag_in[wdm[i]]+rand[*,1]*err_mag_in[wdm[i]]
                        if keyword_set(extrap) then begin
                            log_l_mc = interp2d(log_l_grid_extrap, clr_grid, mag_grid, clr_mc, mag_mc, /extrapolate, /regular)
                            wgood = where(log_l_mc gt -2.0 and log_l_mc lt 7.0, cgood)
                            if cgood ge 2 then fehphot[wdm[i]].err_logl[ii] = stddev(log_l_mc[wgood]) else fehphot[wdm[i]].err_logl[ii] = -999
                        endif else begin
                            log_l_mc = interp2d(log_l_grid, clr_grid, mag_grid, clr_mc, mag_mc, /regular, missing=-999)
                            wgood = where(log_l_mc gt -2.0 and log_l_mc lt 7.0, cgood)
                            if cgood ge 2 then fehphot[wdm[i]].err_logl[ii] = stddev(log_l_mc[wgood]) else fehphot[wdm[i]].err_logl[ii] = -999
                        endelse
                    endelse
                    print, string(27B) + '[1A' + string(double(100*(t+1)) / double(ttot), format='(D5.1)') + '%'
                    t++
                endfor
            endif
            mag_grid -= dm_uniq[jj]
        endfor
    endfor
    wsort = sort(fehphot[0].ages)
    for i=0L,n-1 do begin
        fehphot[i].ages = fehphot[i].ages[wsort]
        fehphot[i].feh = fehphot[i].feh[wsort]
        fehphot[i].teff = fehphot[i].teff[wsort]
        fehphot[i].logg = fehphot[i].logg[wsort]
        fehphot[i].logl = fehphot[i].logl[wsort]
        if noerror eq 0 then begin
            fehphot[i].err_feh = fehphot[i].err_feh[wsort]
            fehphot[i].err_teff = fehphot[i].err_teff[wsort]
            fehphot[i].err_logg = fehphot[i].err_logg[wsort]
            fehphot[i].err_logl = fehphot[i].err_logl[wsort]
        endif
    endfor

    dm = dmin
    return, fehphot
end
