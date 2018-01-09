function sps_interp, z, age
    common sps_spec, sps, spsz, spsage
    if (size(sps))[1] eq 0 then begin
        f = file_search(getenv('ahome')+'fsps/SSP/SSP_Padova_MILES_Kroupa_*.spec')
        sps = sps_read_spec(f)
    endif

    spsz = sps[uniq(sps.zmet, sort(sps.zmet))].zmet
    spsage = sps[uniq(sps.agegyr, sort(sps.agegyr))].agegyr

    if z lt min(spsz) or z gt max(spsz) then message, 'Metallicity is off grid.'
    if age lt min(spsage) or age gt max(spsage) then message, 'Age is off grid.'

    nspsz = n_elements(spsz)
    nspsage = n_elements(spsage) 

    wz = value_locate(spsz, z)
    if spsz[wz] eq z then begin
        iz = wz
        dz = 1d
        nz = 1
    endif else begin
        if wz eq nspsz-1 then wz -= 1
        iz = [wz, wz+1]
        dz = abs(reverse(spsz[iz])-z)
        nz = 2
    endelse

    wa = value_locate(spsage, age)
    if spsage[wa] eq age then begin
        ia = wa
        da = 1d
        na = 1
    endif else begin
        if wa eq nspsage-1 then wa -= 1
        ia = [wa, wa+1]
        da = abs(reverse(spsage[ia])-age)
        na = 2
    endelse

    spsspec = dblarr(5252)
    for i=0,nz-1 do begin
        for j=0,na-1 do begin
            w = where(sps.zmet eq spsz[iz[i]] and sps.agegyr eq spsage[ia[j]], c)
            if c eq 0 then message, 'Spectrum not found.'
            if i eq 0 and j eq 0 then begin
                lambda = sps[w].lambda
                spsspeci = sps[w].spec
            endif else spsspeci = interpol(sps[w].spec, sps[w].lambda, lambda)
            spsspec += dz[i]*da[j]*spsspeci
        endfor
    endfor
    spsspec /= (nz eq 2 ? (spsz[iz[1]]-spsz[iz[0]]) : 1d) $
                *(na eq 2 ? (spsage[ia[1]]-spsage[ia[0]]) : 1d)

    return, {lambda:lambda, spec:spsspec}
end
