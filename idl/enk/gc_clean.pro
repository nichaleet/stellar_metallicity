function gc_clean, moogify, zero=zero, agb=agb
    gc = get_gc()
    n = n_elements(moogify)
    w = bytarr(n)
    masks = moogify[uniq(moogify.mask, sort(moogify.mask))].mask
    nmasks = n_elements(masks)
    for i=0,nmasks-1 do begin
        wm = where(strtrim(moogify.mask, 2) eq strtrim(masks[i], 2))
        wg = where(strtrim(gc.mask, 2) eq strtrim(masks[i], 2), c)
        if c ne 1 then wg = where(strtrim(gc.mask, 2) eq strtrim(masks[i], 2)+'r', c)
        if c ne 1 then begin
            message, 'I could not find mask '+masks[i]+' in the GC structure.';, /info
            w[wm] = 1
        endif
        readcol, getenv('M31')+'specabund/moogify/cmd/cmd_select_'+strtrim(gc[wg].cluster, 2)+(keyword_set(agb) ? '_agb' : '')+'.dat', color, mag, skipline=1, format='D,D', /silent
        if gc[wg].cluster eq 'N5024' or gc[wg].cluster eq 'N7492' then begin
            w[wm] = keyword_set(zero) ? inside(moogify[wm].b0-moogify[wm].v0, moogify[wm].v0, color, mag) : inside(moogify[wm].b-moogify[wm].v, moogify[wm].v, color, mag)
        endif else begin
            w[wm] = keyword_set(zero) ? inside(moogify[wm].v0-moogify[wm].i0, moogify[wm].v0, color, mag) : inside(moogify[wm].v-moogify[wm].i, moogify[wm].v, color, mag)
        endelse
    endfor
    return, w
end
