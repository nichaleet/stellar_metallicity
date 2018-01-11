function get_harris
    harris = {id:' ', name:' ', cluster:' ', ra:0d, dec:0d, l:0d, b:0d, r_sun:0d, r_gc:0d, x:0d, y:0d, z:0d, feh:0d, wt:0, ebv:0d, v_hb:0d, dm:0d, v_t:0d, m_v:0d, ub:0d, bv:0d, vr:0d, vi:0d, spt:' ', ellip:0d, v_r:0d, v_rerr:0d, v_lsr:0d, sigma_v:0d, sigma_verr:0d, c:0d, cflag:' ', r_c:0d, r_h:0d, mu_v:0d, rho_0:0d, logtc:0d, logth:0d}
    harris = replicate(harris, 500)
    n = 0L

    id = ' '
    name = ' '
    radec = ' '
    spt = ' '
    cflag = ' '
    wt = 0

    openr, lun1, getenv('M31')+'gc/harris1.txt', /get_lun
    skip_lun, lun1, 16, /lines
    while ~eof(lun1) do begin
        readf, lun1, id, name, radec, l, b, r_sun, r_gc, x, y, z, format='(A11,1X,A12,1X,A25,2(2X,D6),1X,5(1X,D5))'
        harris[n].id = strtrim(id, 2)
        harris[n].name = strtrim(name, 2)
        get_coords, coords, instring=radec
        harris[n].ra = coords[0]*15.
        harris[n].dec = coords[1]
        harris[n].l = l
        harris[n].b = b
        harris[n].r_sun = r_sun
        harris[n].r_gc = r_gc
        harris[n].x = x
        harris[n].y = y
        harris[n].z = z
        n++
    endwhile
    close, lun1
    free_lun, lun1
    harris = harris[0:n-1]

    openr, lun2, getenv('M31')+'gc/harris2.txt', /get_lun
    skip_lun, lun2, 19, /lines
    while ~eof(lun2) do begin
        readf, lun2, id, feh, wt, ebv, v_hb, dm, v_t, m_v, ub, bv, vr, vi, spt, ellip, format='(A11,2X,D5,1X,I2,1X,4(1X,D5),1X,D6,1X,4(1X,D5),2X,A4,2X,D4)'
        w = where(strtrim(harris.id, 2) eq strtrim(id, 2))
        harris[w].feh = feh
        harris[w].wt = wt
        harris[w].ebv = ebv
        harris[w].v_hb = v_hb
        harris[w].dm = dm
        harris[w].v_t = v_t
        harris[w].m_v = m_v
        harris[w].ub = ub
        harris[w].bv = bv
        harris[w].vr = vr
        harris[w].vi = vi
        harris[w].spt = spt
        harris[w].ellip = ellip
    endwhile
    close, lun2
    free_lun, lun2

    openr, lun3, getenv('M31')+'gc/harris3.txt', /get_lun
    skip_lun, lun3, 20, /lines
    while ~eof(lun3) do begin
        readf, lun3, id, v_r, v_rerr, v_lsr, sigma_v, sigma_verr, c, cflag, r_c, r_h, mu_v, rho_0, logtc, logth, format='(A11,2X,D6,2X,D4,2X,D6,4X,D4,2X,D4,3X,D4,2X,A2,2(2X,D4),3X,D5,2(2X,D5),1X,D5)'
        w = where(strtrim(harris.id, 2) eq strtrim(id, 2))
        harris[w].v_r = v_r
        harris[w].v_rerr = v_rerr
        harris[w].v_lsr = v_lsr
        harris[w].sigma_v = sigma_v
        harris[w].sigma_verr = sigma_verr
        harris[w].c = c
        harris[w].cflag = cflag
        harris[w].r_c = r_c
        harris[w].r_h = r_h
        harris[w].mu_v = mu_v
        harris[w].rho_0 = rho_0
        harris[w].logtc = logtc
        harris[w].logth = logth
    endwhile
    close, lun3
    free_lun, lun3

    for i=0,n-1 do begin
        cluster = strcompress(harris[i].id, /remove_all)
        if strmid(cluster, 0, 3) eq 'NGC' then cluster = 'N'+strmid(cluster, 3)
        harris[i].cluster = cluster
        if strmid(harris[i].name, 0, 2) eq 'M ' then harris[i].name = 'M'+strmid(harris[i].name, 2)
        if strtrim(harris[i].name, 2) eq 'AvdB' then harris[i].name = harris[i].id
        if strtrim(harris[i].name, 2) eq '' then harris[i].name = harris[i].id
    endfor

    return, harris
end
