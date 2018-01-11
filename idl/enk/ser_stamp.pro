;+
; NAME:
;    ser_stamp
;
; PURPOSE:
;    to create an ds9 plot with postage stamps of a given list of objects
;      - note: all objects have to be in the same pointing (ie. field '12')
;
; CATEGORY:
;    photometric analysis
;
;
; CALLING SEQUENCE:
;	enk_stamp, struct,  ds9, stampsz = stampsz, $
;		outfile = outfile, mask=mask
;
; INPUTS:
;    struct = input structure containing list of objects you want postage
;                    stamps of - ex: photo or targ
; 
; OPTIONAL INPUTS:
;    none
;	
; KEYWORD PARAMETERS:
;    stampsz = sets the pixel width of each postage stamp - default 50 pix
;    outfile = string name of fits file to save the composite large image to
;    mask = Maskname, string, used to determine PA of slit
;
; OUTPUTS:
;    creates an atvplot of one large image comprised of all the
;            postage stamps, with optional ellipses drawn for each object
;
;
; MODIFICATION HISTORY:
;    enk 28jul05 - fixed slit length and position
;    npk 15jan05 - Significant rewrite to work with ds9 + provide slitimage
;    alc 28feb02
;    alc 07mar02 - loop through quilts only once
;    alc 22mar02 - fixed sorting of quilts - now works for "targets.fits"
;    alc 08dec02 - added the color and allcolor keywords
;-

pro drawslit, cat, struct, ds9, cfh=cfh, acs=acs
    ;catch, error_status
    error_status = 0
    if error_status ne 0 then begin
        message, 'Could not draw slit for some reason!', /info
        return
    endif
    mn = strcompress(cat.mask,/rem)
    ip = mrdfits(getenv('D2_RESULTS')+'/../1HSmasks/' + mn + '/1HSmask.' + mn + '.fits', 1)
    mm = mrdfits(getenv('D2_RESULTS')+'/../1HSmasks/' + mn + '/1HSmask.' + mn + '.fits', 2)

    r = where(mm.objno eq struct.objno, count)

    if count eq 0 then begin
        message, 'Could not find ' + string(struct.objno) + ' in photometry file.  Can not overplot slits.',/info
        return
    endif

    ra_obj = mm[r].ra & dec_obj = mm[r].dec
    o = mm[r]

    dx = o.slitwidth/3600.
    dy = (o.topdist+o.botdist)/3600.
    pa = o.slitpa+ip.pa
    angle = string(pa,format='(F15.9)')
    pa = pa * !DTOR
    dec_obj_rad = dec_obj * !DTOR

    ra = ra_obj + 0.5*(o.topdist-o.botdist)/3600.*sin(pa)/cos(dec_obj_rad)
    dec = dec_obj + 0.5*(o.topdist-o.botdist)/3600.*cos(pa)
    sn = string(o.slitn)

    ;x1 = ra-0.5*dy*sin(pa)/cos(dec_obj_rad)
    ;y1 = dec-0.5*dy*cos(pa)
    ;x2 = ra+0.5*dy*sin(pa)/cos(dec_obj_rad)
    ;y2 = dec+0.5*dy*cos(pa)
    ;slope_a = (y2-y1)/(x2-x1)
    ;x1 = string(x1,format='(F15.9)')
    ;y1 = string(y1,format='(F15.9)')
    ;x2 = string(x2,format='(F15.9)')
    ;y2 = string(y2,format='(F15.9)')
    ;print, x1 + y1 + x2 + y2
    ;print, 'slope_a: ' + string(slope_a,format='(F15.9)')
    ;ds9->pipecmd, 'fk5; line' + x1 + y1 + x2 + y2 + '# color=blue', 'regions'

    ;x3 = ra+0.5*dx*cos(pa)/cos(dec_obj_rad)
    ;y3 = dec-0.5*dx*sin(pa)
    ;x4 = ra-0.5*dx*cos(pa)/cos(dec_obj_rad)
    ;y4 = dec+0.5*dx*sin(pa)
    ;slope_b = (y4-y3)/(x4-x3)
    ;x3 = string(x3,format='(F15.9)')
    ;y3 = string(y3,format='(F15.9)')
    ;x4 = string(x4,format='(F15.9)')
    ;y4 = string(y4,format='(F15.9)')
    ;print, x3 + y3 + x4 + y4
    ;print, 'slope_b: ' + string(slope_b,format='(F15.9)')
    ;print, 'slope_a * slope_b: ' + string(slope_a*slope_b,format='(F15.9)')
    ;ds9->pipecmd, 'fk5; line' + x3 + y3 + x4 + y4 + '# color=blue', 'regions'

    ra = string(ra,format='(F15.9)')
    dec = string(dec,format='(F15.9)')        
    dx = string(dx,format='(F15.9)')
    dy = string(dy,format='(F15.9)')
    pa = string(pa,format='(F15.9)')

    ;print, 'angle: ' + angle
    ;print, 'pa: ' + pa

    ds9->pipecmd, 'fk5; box' + ra + dec + dx + dy + angle + '# color=red', 'regions'
    ;ds9->pipecmd, 'fk5; box' + ra + dec + dx + dy + angle + '# text={' + sn + '} color=red', 'regions'
    ;ds9->pipecmd, 'fk5; box' + ra + dec + dx + dy + string(29) + '# color=blue', 'regions'

    ra_obj = string(ra_obj,format='(F15.9)')
    dec_obj = string(dec_obj,format='(F15.9)')
    ds9->pipecmd, 'fk5; point' + ra_obj + dec_obj + '# point=cross color=blue', 'regions'
    ;ds9->pipecmd, 'fk5; point' + ra + dec + '# point=cross color=green', 'regions'

    yra = string(cat.serra,format='(F15.9)')
    ydec = string(cat.serdec,format='(F15.9)')
    ;print, yra, ydec
    ds9->pipecmd, 'fk5; point' + yra + ydec + '# point=cross color=green', 'regions'


    wcfh = where(cat.cfh_type eq 2, cwcfh)
    if keyword_set(cfh) and cwcfh gt 0 then begin
        for i=0,cwcfh-1 do begin
            ;ytoradec, mn, struct.objno, yoffset, yra, ydec
            yra = string(cat.cfh_ra[wcfh[i]],format='(F15.9)')
            ydec = string(cat.cfh_dec[wcfh[i]],format='(F15.9)')
            ;print, yra, ydec
            ds9->pipecmd, 'fk5; point' + yra + ydec + '# point=circle color=yellow', 'regions'
        endfor
    endif
    wcfh = where(cat.cfh_type eq 3, cwcfh)
    if keyword_set(cfh) and cwcfh gt 0 then begin
        for i=0,cwcfh-1 do begin
            ;ytoradec, mn, struct.objno, yoffset, yra, ydec
            yra = string(cat.cfh_ra[wcfh[i]],format='(F15.9)')
            ydec = string(cat.cfh_dec[wcfh[i]],format='(F15.9)')
            ;print, yra, ydec
            ds9->pipecmd, 'fk5; point' + yra + ydec + '# point=circle color=green', 'regions'
        endfor
    endif
    wacs = where(cat.acs_tile ne 0 and cat.acs_id ne 0, cwacs)
    if keyword_set(acs) and cwacs gt 0 then begin
        for i=0,cwacs-1 do begin
            ;ytoradec, mn, struct.objno, yoffset, yra, ydec
            yra = string(cat.acs_ra[wacs[i]],format='(F15.9)')
            ydec = string(cat.acs_dec[wacs[i]],format='(F15.9)')
            ;print, yra, ydec
            ds9->pipecmd, 'fk5; point' + yra + ydec + '# point=circle color=green', 'regions'
        endfor
    endif

    ;ds9->cmd, 'regions save ' + mn + '.' + strcompress(sn,/rem) + '.reg'

END

pro update_astrometry, a, d, ast, ds9
    if (size(ast))[0] eq 0 then begin
        return
    endif
    ast.crpix = [0,0]
    ast.crval = [a,d]
    hdr = strarr(1)
    putast, hdr, ast
    s = strcompress(string(randomu(s,1)),/rem)
    openw, 1, '/tmp/' + s[0] + 'wcs.tmp'
    for i = 0, n_elements(hdr)-1 do begin
        printf, 1, hdr[i]
    endfor
    close,1
    ds9->cmd, 'wcs append file /tmp/' + s[0] + 'wcs.tmp'
    file_delete, '/tmp/' + s[0] + 'wcs.tmp'
END

;pro enk_stamp, struct,  ds9, stampsz = stampsz, $
;   outfile = outfile, mask=mask, yoffset=yoffset

pro ser_stamp, cat, ds9, stampsz=stampsz, outfile=outfile
    objnos = strcompress(cat.tgtobjno, /rem)
    field = strmid(objnos,0,2)
    pp = getphot(field)
    struct = pp[where(pp.objno eq cat.tgtobjno)]

    ; pixel width of each postage_stamp
    if NOT keyword_set(stampsz) then stampsz = 200

    ; find which quilt it is on
    quilt = quilt_number(struct.xs, field, x, /string)

    image = intarr(3,stampsz, stampsz)

    ; read in R band quilt image
    field1 = strmid(field, 0, 1)
    field2 = strmid(field, 1, 1)
    path = 'field'+field1+'/'+field2
    directory=getenv('IDLDEEPPHOTOMETRY_DIR')  + path  + '/quilt/'        

    rimage=mrdfits(directory + 'R.'+quilt+'.fits.Z',/silent, 0, header)
    extast, header, rast

    bimage=mrdfits(directory + 'B.'+quilt+'.fits.Z',/silent, 0, header)
    extast, header, bast

    iimage=mrdfits(directory + 'I.'+quilt+'.fits.Z',/silent, 0, header)
    extast, header, iast

    ; get postage stamps 
    sz = stampsz
    qltsz = 1024

    ;ds9->cmd, 'frame delete all'
    objname = strcompress(struct.objno, /rem)

    tmp = quilt_number(struct.xs, field, x, /string)
    x[0] = round(x[0]) & x[1] = round(x[1])
    xp = x-sz/2-1 ; xp is in the coordinates of the astrometry

    ; put quilt in a larger image to take care of cases where stamp 
    ; hits edge of quilt

    largeim = fltarr(qltsz+sz,qltsz+sz)
    largeim[sz/2-1:sz/2-1+qltsz-1, sz/2-1:sz/2-1+qltsz-1] = bimage
    
    ; extract the B stamp
    ds9->cmd, 'frame 0'
    stamp = largeim[x[0]-1:x[0]+sz-2, x[1]-1:x[1]+sz-2]
    xstart = 0 & ystart = 0
    ds9->image,stamp
    ds9->cmd, 'scale mode zscale'
    ds9->cmd, 'zoom to 4'
    ds9->cmd, 'wcs fk5'        
    ds9->cmd, 'wcs align yes'
    ds9->cmd, 'wcs format degrees'
    xy2ad, xp[0]-1, xp[1]-1, bast, a, d
    update_astrometry, a, d, bast, ds9

    largeim[sz/2-1:sz/2-1+qltsz-1, sz/2-1:sz/2-1+qltsz-1] = rimage
    ; extract the R stamp
    ds9->cmd, 'frame 1'
    stamp = largeim[x[0]-1:x[0]+sz-2, x[1]-1:x[1]+sz-2]
    xstart = 0 & ystart = 0
    ds9->image,stamp
    ds9->cmd, 'scale mode zscale'
    xy2ad, xp[0]-1, xp[1]-1, rast, a, d
    update_astrometry, a, d, rast, ds9

    drawslit, cat, struct, ds9, /cfh

    largeim[sz/2-1:sz/2-1+qltsz-1, sz/2-1:sz/2-1+qltsz-1] = iimage
    ; extract the I stamp
    ds9->cmd, 'frame 2'
    stamp = largeim[x[0]-1:x[0]+sz-2, x[1]-1:x[1]+sz-2]
    xstart = 0
    ystart = 0
    ds9->image,stamp
    ds9->cmd, 'scale mode zscale'
    xy2ad, xp[0]-1, xp[1]-1, iast, a, d
    update_astrometry, a, d, iast, ds9                

    ;draw ellipses on the I stamp
    xoffset = (sz-1)/2.
    yoffset = (sz-1)/2.
    a = 3*struct.rg
    e2 = struct.e2
    b = (1-e2)*a
    ;ds9->ellipse, xoffset, yoffset, theta=struct.pa, b,a,
    ;color="blue"
    ras = string(struct.ra,format='(F15.9)')
    decs = string(struct.dec,format='(F15.9)')
    as = string(a*0.207/3600.,format='(F15.9)')
    bs = string(b*0.207/3600.,format='(F15.9)')
    pas = string(struct.pa,format='(F15.9)')
    ds9->pipecmd, 'fk5; ellipse' + ras + decs + bs + as + pas + '# color=blue', 'regions'
    w = where(cat.cfh_type eq 2, w2)
    for i=0,w2-1 do begin
        objno2 = strcompress(cat.cfh_objno[w[i]], /rem)
        field2 = strmid(objno2,0,2)
        pp = getphot(field2)
        struct2 = pp[where(pp.objno eq cat.cfh_objno[w[i]])]
        a = 3*struct2.rg
        e2 = struct2.e2
        b = (1-e2)*a
        ;ds9->ellipse, xoffset, yoffset, theta=struct2.pa, b,a, color="yellow"
        ras = string(struct2.ra,format='(F15.9)')
        decs = string(struct2.dec,format='(F15.9)')
        as = string(a*0.207/3600.,format='(F15.9)')
        bs = string(b*0.207/3600.,format='(F15.9)')
        pas = string(struct2.pa,format='(F15.9)')
        ds9->pipecmd, 'fk5; ellipse' + ras + decs + bs + as + pas + '# color=yellow', 'regions'
    endfor
    w = where(cat.cfh_type eq 3, w3)
    for i=0,w3-1 do begin
        objno3 = strcompress(cat.cfh_objno[w[i]], /rem)
        field3 = strmid(objno3,0,2)
        pp = getphot(field3)
        struct3 = pp[where(pp.objno eq cat.cfh_objno[w[i]])]
        a = 3*struct3.rg
        e2 = struct3.e2
        b = (1-e2)*a
        ;ds9->ellipse, xoffset, yoffset, theta=struct3.pa, b,a, color="green"
        ras = string(struct3.ra,format='(F15.9)')
        decs = string(struct3.dec,format='(F15.9)')
        as = string(a*0.207/3600.,format='(F15.9)')
        bs = string(b*0.207/3600.,format='(F15.9)')
        pas = string(struct3.pa,format='(F15.9)')
        ds9->pipecmd, 'fk5; ellipse' + ras + decs + bs + as + pas + '# color=green', 'regions'
    endfor


    vfilename = '/net/marc/raid2/egs_stamps/deep2/'+strmid(objname,0,5)+'/'+objname+'.acs.v.fits.gz'
    ;print, vfilename
    if file_test(vfilename) then begin
        ds9->cmd, 'frame 3'
        stamp = mrdfits(vfilename, 0, hdr)
        xstart = 0
        ystart = 0
        ds9->image,stamp
        ds9->update_wcs, hdr
        ds9->cmd, 'scale mode zscale'
    endif

    ifilename = '/net/marc/raid2/egs_stamps/deep2/'+strmid(objname,0,5)+'/'+objname+'.acs.i.fits.gz'
    if file_test(ifilename) then begin
        ds9->cmd, 'frame 4'
        stamp = mrdfits(ifilename, 0)
        xstart = 0
        ystart = 0
        ds9->image,stamp
        ds9->update_wcs, hdr
        ds9->cmd, 'scale mode zscale'
        drawslit, cat, struct, ds9, /acs
    endif
    ds9->cmd, 'frame first'
    ds9->cmd, 'match frames wcs'

    if keyword_set(outfile) then writefits, outfile, image
end

