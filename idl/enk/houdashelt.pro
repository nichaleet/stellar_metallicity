function houdashelt, color_in, logg_in, feh_in, ub=ub, bv=bv, vr=vr, vi=vi, vj=vj, vk=vk, jk=jk, hk=hk, bc=bc
    grid = {teff:0d, logg:0d, feh:0d, ub:0d, bv:0d, vr:0d, vi:0d, vj:0d, vk:0d, jk:0d, hk:0d, bck:0d}
    grid = replicate(grid, 10000)

    i = 0L
    ;openr, lun, getenv('CALTECH')+'houdashelt/hou00a_tab5.dat', /get_lun
    ;while ~eof(lun) do begin
    ;    readf, lun, teff_read, logg_read, feh_read, ub_read, bv_read, vr_read, vi_read, vk_read, jk_read, hk_read, bck_read, format='(I4,3X,F4.2,2X,F5.2,11X,7(F5.3,3X),32X,F5.3)'
    ;    if bck_read le 0 then continue
    ;    grid[i].teff = teff_read
    ;    grid[i].logg = logg_read
    ;    grid[i].feh = feh_read
    ;    grid[i].ub = ub_read
    ;    grid[i].bv = bv_read
    ;    grid[i].vr = vr_read
    ;    grid[i].vi = vi_read
    ;    grid[i].vk = vk_read
    ;    grid[i].jk = jk_read
    ;    grid[i].hk = hk_read
    ;    grid[i].bck = bck_read
    ;    i++
    ;endwhile
    ;close, lun
    ;free_lun, lun

    openr, lun, getenv('CALTECH')+'houdashelt/hou00b_tab5.dat', /get_lun
    while ~eof(lun) do begin
        readf, lun, teff_read, logg_read, feh_read, ub_read, bv_read, vr_read, vi_read, vk_read, jk_read, hk_read, bck_read, format='(I4,3X,F4.2,3X,F5.2,2X,4(F5.2,11X),3(F5.3,3X),32X,F5.3)'
        if bck_read le 0 then continue
        grid[i].teff = teff_read
        grid[i].logg = logg_read
        grid[i].feh = feh_read
        grid[i].ub = ub_read
        grid[i].bv = bv_read
        grid[i].vr = vr_read
        grid[i].vi = vi_read
        grid[i].vj = vk_read - jk_read
        grid[i].vk = vk_read
        grid[i].jk = jk_read
        grid[i].hk = hk_read
        grid[i].bck = bck_read
        i++
    endwhile
    close, lun
    free_lun, lun
    grid = grid[0:i-1]

    if keyword_set(ub) then ic = 3
    if keyword_set(bv) then ic = 4
    if keyword_set(vr) then ic = 5
    if keyword_set(vi) then ic = 6
    if keyword_set(vj) then ic = 7
    if keyword_set(vk) then ic = 8
    if keyword_set(jk) then ic = 9
    if keyword_set(hk) then ic = 10

    ulogg = grid[uniq(grid.logg, sort(grid.logg))].logg
    ufeh = grid[uniq(grid.feh, sort(grid.feh))].feh

    wlogg = value_locate(ulogg, logg_in)
    if wlogg lt 0 then begin
        logg_in = min(grid.logg)
        wlogg = 0
    endif
    if ulogg[wlogg] eq logg_in then begin
        ilogg = wlogg
        dlogg = 1d
        ng = 1
    endif else begin
        if wlogg eq n_elements(ulogg)-1 then wlogg -= 1
        ilogg = [wlogg, wlogg+1]
        dlogg = abs(reverse(ulogg[ilogg])-logg_in)
        ng = 2
    endelse

    wfeh = value_locate(ufeh, feh_in)
    if ufeh[wfeh] eq feh_in then begin
        ifeh = wfeh
        dfeh = 1d
        nf = 1
    endif else begin
        if wfeh eq n_elements(ufeh)-1 then wfeh -= 1
        ifeh = [wfeh, wfeh+1]
        dfeh = abs(reverse(ufeh[ifeh])-feh_in)
        nf = 2
    endelse
    
    for j=0,ng-1 do begin
        for m=0,nf-1 do begin
            w = where(grid.logg eq ulogg[ilogg[j]] and grid.feh eq ufeh[ifeh[m]])
            teffi = grid[w].teff
            colori = grid[w].(ic)
            bcki = grid[w].bck
            if (size(teff_interp))[1] eq 0 then begin
                teff_interp = dlogg[j]*dfeh[m]*teffi
                color_interp = dlogg[j]*dfeh[m]*colori
                bck_interp = dlogg[j]*dfeh[m]*bcki
            endif else begin
                teff_interp += dlogg[j]*dfeh[m]*teffi
                color_interp += dlogg[j]*dfeh[m]*colori
                bck_interp += dlogg[j]*dfeh[m]*bcki
            endelse
        endfor
    endfor
    teff_interp /= (ng eq 2 ? (ulogg[ilogg[1]]-ulogg[ilogg[0]]) : 1d) $
                    *(nf eq 2 ? (ufeh[ifeh[1]]-ufeh[ifeh[0]]) : 1d)
    color_interp /= (ng eq 2 ? (ulogg[ilogg[1]]-ulogg[ilogg[0]]) : 1d) $
                    *(nf eq 2 ? (ufeh[ifeh[1]]-ufeh[ifeh[0]]) : 1d)
    bck_interp /= (ng eq 2 ? (ulogg[ilogg[1]]-ulogg[ilogg[0]]) : 1d) $
                  *(nf eq 2 ? (ufeh[ifeh[1]]-ufeh[ifeh[0]]) : 1d)

    teff = interpol(teff_interp, color_interp, color_in)
    bck = interpol(bck_interp, color_interp, color_in)
    bc = keyword_set(vk) ? bck - color_in : bck
    return, teff
end
