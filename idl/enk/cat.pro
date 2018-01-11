function cat, espec, objnames, v, verr, fit=fit, dm=dm, errdm=errdm
    if ~keyword_set(fit) then fit = 2
    if ~keyword_set(dm) then message, 'You must specify the distance modulus (keyword DM).'
    if ~keyword_set(errdm) then errdm = 0.0

    nobj = n_elements(objnames)
    cat = {objname:' ', mv:0d, mverr:0d, a:dblarr(3, 5), ewcat:dblarr(3), ewcaterr:dblarr(3), fehcat:0d, fehcaterr:0d}
    cat = replicate(cat, nobj)
    cat.objname = objnames

    wcat1 = where(espec.atomic eq 20 and espec.ion eq 1 and espec.fit eq fit and round(espec.lambda) eq 8498)
    wcat2 = where(espec.atomic eq 20 and espec.ion eq 1 and espec.fit eq fit and round(espec.lambda) eq 8542)
    wcat3 = where(espec.atomic eq 20 and espec.ion eq 1 and espec.fit eq fit and round(espec.lambda) eq 8662)
    match, strtrim(cat.objname, 2), strtrim(espec[wcat1].objname, 2), w1, w2
    cat[w1].ewcat[0] = espec[wcat1[w2]].ew*0.75
    cat[w1].ewcaterr[0] = espec[wcat1[w2]].ewerr*0.75
    for i=0,n_elements(w1)-1 do cat[w1[i]].a[0,*] = espec[wcat1[w2[i]]].a
    match, cat.objname, espec[wcat2].objname, w1, w2
    cat[w1].ewcat[1] = espec[wcat2[w2]].ew*0.75
    cat[w1].ewcaterr[1] = espec[wcat2[w2]].ewerr*0.75
    for i=0,n_elements(w1)-1 do cat[w1[i]].a[1,*] = espec[wcat2[w2[i]]].a
    match, cat.objname, espec[wcat3].objname, w1, w2
    cat[w1].ewcat[2] = espec[wcat3[w2]].ew*0.75
    cat[w1].ewcaterr[2] = espec[wcat3[w2]].ewerr*0.75
    for i=0,n_elements(w1)-1 do cat[w1[i]].a[2,*] = espec[wcat3[w2[i]]].a
    for i=0,n_elements(cat)-1 do begin
        cat[i].mv = v[i] - dm
        cat[i].mverr = sqrt(verr[i]^2. + errdm^2.)
        if cat[i].ewcat[1] gt 0.0 and cat[i].ewcaterr[1] gt 0.0 and cat[i].ewcat[2] gt 0.0 and cat[i].ewcaterr[2] gt 0.0 then begin
            ewsum = 1d-3*(cat[i].ewcat[1]+cat[i].ewcat[2])
            ewsumerr = 1d-3*sqrt((cat[i].ewcaterr[1])^2.+(cat[i].ewcaterr[2])^2.)
            cat[i].fehcat = -2.90 + 0.187*cat[i].mv + 0.422*ewsum - 0.882*ewsum^(-1.5) + 0.0133*ewsum*cat[i].mv
            cat[i].fehcaterr = sqrt(((0.187 + 0.0133*ewsum)*cat[i].mverr)^2. + ((0.422 + 1.5*0.822*ewsum^(-2.5) + 0.0133*cat[i].mv)*ewsumerr)^2.)
        endif else begin
            cat[i].fehcat = -999.0
            cat[i].fehcaterr = -999.0
        endelse
    endfor
    return, cat
end
